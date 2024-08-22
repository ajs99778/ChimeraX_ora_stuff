import os
import subprocess

from chimerax.atomic import AtomicStructure

import numpy as np

from AaronTools.const import ELEMENTS

from ora_stuff.mouse_modes import (
    SetHalfBond,
    SetSingleBond,
    SetAromaticBond,
    SetPartialDoubleBond,
    SetDoubleBond,
    SetTripleBond,
)

bo_to_mol_map = {
    0.5: 8,
    1.0: 1,
    "ar": 4,
    1.5: 5,
    2.0: 2,
    3.0: 3,
}

mol_to_bo_map = {
    8: "half",
    1: "single",
    5: "partial double",
    4: "aromatic",
    2: "double",
    3: "triple",
}

def open_sdf(session, stream, filename, **kw):
    try:
        all_coordsets = []
        elements = []
        all_bondsets = []
        line = stream.readline()
        n_atoms = 0
        n_bonds = 0
        while line:
            if "v2000" in line.lower():
                # read according to the MDL Informatics Systems prescription
                # from 2003
                this_coordset = []
                this_bondset = []
                n_atoms = int(line[:3])
                n_bonds = int(line[3:6])
                for i in range(0, n_atoms):
                    # columns 1-10 define x position, and there should be
                    # 4 digits after the decimal
                    # columns 11-20 are for the y position
                    # columns 21-30 are for the z position
                    # column 31 is empty
                    # columns 32-34 are for the element symbol
                    line = stream.readline()
                    coords = [float(line[0:11]), float(line[11:21]), float(line[21:31])]
                    this_coordset.append(coords)
                    if all_coordsets:
                        continue
                    ele = line[31:34].strip()
                    elements.append(ele)
                
                for i in range(0, n_bonds):
                    line = stream.readline()
                    a1 = int(line[:3])
                    a2 = int(line[3:6])
                    order = int(line[6:9])
                    this_bondset.append([a1 - 1, a2 - 1, order])
                
                all_coordsets.append(this_coordset)
                all_bondsets.append(this_bondset)
            
            elif "v3000" in line.lower():
                this_coordset = []
                this_bondset = []
                while "M END" not in line and line.strip():
                    line = stream.readline()
                    if "COUNTS" in line:
                        n_atoms = int(line.split()[3])
                        n_bonds = int(line.split()[4])

                    if "BEGIN ATOM" in line:
                        for i in range(0, n_atoms):
                            line = stream.readline()
                            x, y, z = [float(x) for x in line.split()[4:7]]
                            this_coordset.append([x, y, z])
                            if all_coordsets:
                                continue
                            ele = line.split()[3].strip()
                            elements.append(ele)
                    
                    if "BEGIN BOND" in line:
                        for i in range(0, n_bonds):
                            line = stream.readline()
                            info = line.split()
                            a1 = int(info[4])
                            a2 = int(info[5])
                            order = int(info[3])
                            this_bondset.append([a1 - 1, a2 - 1, order])
                
                all_coordsets.append(this_coordset)
                all_bondsets.append(this_bondset)
            
            else:
                line = stream.readline()

        stream.close()

        struc = AtomicStructure(session)
        struc.name = filename
        res = struc.new_residue("UNK", "a", 1)
        ele_counts = dict()
        for i, ele in enumerate(elements):
            try:
                ele = ELEMENTS[int(ele)]
            except ValueError:
                pass
            except IndexError:
                raise RuntimeError(
                    "error while trying to convert atomic number to symbol: %s" % ele
                )
            ele_counts.setdefault(ele, 0)
            ele_counts[ele] += 1
            name = "%s%i" % (ele, ele_counts[ele])
            atom = struc.new_atom(name, ele)
            atom.coord = np.array(all_coordsets[0][i])
            atom.serial_number = i + 1
            res.add_atom(atom)
        
        struc.add_coordsets(np.array(all_coordsets), replace=True)
        for cs_id, bondset in zip(struc.coordset_ids, all_bondsets):
            for bond in bondset:
                a1, a2, order = bond
                try:
                    pbg = struc.pseudobond_group(
                        mol_to_bo_map[order],
                        create_type=2
                    )
                    pbg.new_pseudobond(struc.atoms[a1], struc.atoms[a2], cs_id)
                except KeyError:
                    session.logger.warning("unknown bond order %i" % order)
                    continue
        
        for bo_pbg, style in zip(
            ["half", "single", "aromatic", "partial double", "double", "triple"],
            [SetHalfBond, SetSingleBond, SetAromaticBond, SetPartialDoubleBond, SetDoubleBond, SetTripleBond]
        ):
            pbg = struc.pseudobond_group(
                bo_pbg,
                create_type=None,
            )
            if pbg is None:
                continue
            pbg.dashes = style.dashes
            pbg.dashes += 2
            pbg.dashes -= 2
            pbg.color = style.color
        
        struc.active_coordset_id = struc.coordset_ids[0]
        
        return [struc], "opened file"
        
    except Exception as e:
        stream.close()
        raise e

def save_sdf(
    session,
    path,
    model=None,
    style="V3000",
    coordsets=True,
):

    if coordsets:
        write_coordsets = model.coordset_ids
    else:
        write_coordsets = [model.active_coordset_id]
    
    # write the mol3 file as we go through each set of coordinates
    with open(path, "w") as f:
        for cs_id in write_coordsets:
            f.write("%s\n" % model.name)
            f.write("ORA stuff for ChimeraX\n")
            f.write("\n")

            coords = model.coordset(cs_id).xyzs
            # squared distance for each pair of atoms
            
            this_bonds = []
            ndx = {a: i for i, a in enumerate(model.atoms)}
            for pbg_name, order in zip(
                ["half", "single", "aromatic", "partial double", "double", "triple"],
                [0.5, 1.0, "ar", 1.5, 2.0, 3.0]
            ):
                pbg = model.pseudobond_group(
                    pbg_name,
                    create_type=None
                )
                if not pbg:
                    continue
                
                for bond in pbg.get_pseudobonds(cs_id):
                    this_bonds.append([ndx[bond.atoms[0]], ndx[bond.atoms[1]], order])

            if style == "V3000":
                f.write("  0  0  0  0  0  0  0  0  0  0  0 V3000\n")
                f.write("M  V30 BEGIN CTAB\n")
                
                f.write("M  V30 COUNTS %i %i 0 0 0\n" % (model.num_atoms, len(this_bonds)))
                f.write("M  V30 BEGIN ATOM\n")
                for i, a in enumerate(model.atoms):
                    f.write("M  V30 %4i %2s %10.6f %10.6f %10.6f 0\n" % (
                        i + 1, a.element.name, *coords[i]
                    ))
                f.write("M  V30 END ATOM\n")
                f.write("M  V30 BEGIN BOND\n")
                for i, bond in enumerate(this_bonds):
                    f.write("M  V30 %3i %i %i %i\n" % (
                        i + 1,
                        bo_to_mol_map[bond[2]],
                        bond[0] + 1,
                        bond[1] + 1
                    ))
                f.write("M  V30 END BOND\n")
                f.write("M  V30 END CTAB\n")
            
            elif style == "V2000":
                f.write("%3i%3i  0  0  0  0  0  0  0  0  0 V2000\n" % (model.num_atoms, len(this_bonds)))
                for i, a in enumerate(model.atoms):
                    f.write("%10.4f%10.4f%10.4f %3s 0%3i  0  0  0  0  0  0  0  0\n" % (
                        *coords[i],
                        a.element.name,
                        0
                    ))
                for bond in this_bonds:
                    f.write("%3i%3i%3i  0  0  0  0\n" % (
                        bond[0] + 1,
                        bond[1] + 1,
                        bo_to_mol_map[bond[2]],
                    ))
                
            else:
                raise NotImplementedError(style)
            
            f.write("M END\n")
            f.write("$$$$\n")

def save_fbx(session, path, model=None, blenderPath=None, scriptOnly=False):
    print(path, model, blenderPath, scriptOnly)
    
    y_cor = 2.2
    
    out = []
    out.append("import os")
    out.append("import bpy")
    out.append("import bmesh")
    out.append("bpy.context.scene.view_settings.view_transform = 'Raw'")
    pad = 1 + int(max(0, np.log10(model.num_atoms)))
    # the spheres for the atoms get a unique label that is the
    # element + atom number (zero padded)
    label_format = "%%s%%0%ii" % pad
    created_materials = []
    for i, a in enumerate(model.atoms):
        ele = a.element.name
        atom_label = label_format % (ele, i + 1)
        # create a sphere for this atom
        out.append(
            "bpy.ops.mesh.primitive_uv_sphere_add(radius=%f, location=(%f,%f,%f))" % (
                a.display_radius, *a.coord
            )
        )
        out.append("ob = bpy.data.objects[\"Sphere\"]")
        out.append("ob.name = " + "\"%s\"" % atom_label)
        if ele not in created_materials:
            out.append("mat_%s = bpy.data.materials.new(name=\"%s\")" % (ele, ele))
            out.append("mat_%s.use_nodes = False" % ele)
            out.append("mat_%s.diffuse_color = (%f,%f,%f,%f)" % (
                ele, (a.color[0] / 255.) ** y_cor, (a.color[1] / 255.) ** y_cor, (a.color[2] / 255.) ** y_cor, a.color[3] / 255.
            ))
            created_materials.append(ele)
            # out.append("bpy.data.materials[\"%s\"].node_tree.nodes[\"Principled BSDF\"].inputs[0].default_value = (%f,%f,%f,%f)" % (
            #     ele, a.color[0] / 255., a.color[1] / 255., a.color[2] / 255., a.color[3] / 255.
            # ))
        # set material and do smooth shading so it looks nice
        out.append("ob.active_material = mat_%s" % ele)
        out.append("bpy.ops.object.shade_smooth()")
    
    # delete the standard cube, light, and camera
    out.append("objs = bpy.data.objects")
    out.append("objs.remove(objs[\"Cube\"], do_unlink=True)")
    out.append("objs.remove(objs[\"Light\"], do_unlink=True)")
    out.append("objs.remove(objs[\"Camera\"], do_unlink=True)")
    
    # fbx file will be exported to ~/chimerax_2_blender/{name}.fbx
    # we make sure that directory exists, but expanduser("~") might
    # have issues on virtual machines or if you're running as admin/root
    # it also might not have issues
    # if it does have issues, they're your issues and not mine
    # figure it out
    out.append("export_path = \"%s\"" % path)
    out.append("os.makedirs(os.path.dirname(export_path), exist_ok=True)")
    out.append("bpy.ops.object.select_all(action='DESELECT')")

    out.append("bpy.ops.export_scene.fbx(filepath=export_path, check_existing=False)")
    # close blender when done
    out.append("import sys; sys.exit(0)")

    # the corresponding python file with all the above code
    # is written to ~/chimerax_2_blender/{name}.py
    script_name, ext = os.path.splitext(path)
    script_name += ".py"
    os.makedirs(os.path.dirname(script_name), exist_ok=True)
    with open(script_name, "w") as f:
        f.write("\n".join(out))
    session.logger.info("saved python file to %s" % script_name)
    
    if scriptOnly:
        return
        
    # run blender on the .py file
    proc = subprocess.Popen(
        [blenderPath, "-P", script_name, "-b"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf8",
    )
    out, err = proc.communicate()
    if out is not None:
        session.logger.info("<pre>" + out + "</pre>", is_html=True)
    if err is not None:
        session.logger.warning("<pre>" + err + "</pre>", is_html=True)