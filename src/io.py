import os
import subprocess

from chimerax.atomic import AtomicStructure

import numpy as np

bo_to_mol_map = {
    0.5: 8,
    1.0: 1,
    1.5: 4,
    2.0: 2,
    3.0: 3,
}

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
                ["half", "single", "aromatic", "double", "triple"],
                [0.5, 1.0, 1.5, 2.0, 3.0]
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