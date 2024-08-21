from chimerax.atomic import PseudobondGroup, AtomsArg, Atoms
from chimerax.core.commands import IntArg, FloatArg, Or, EnumOf, CmdDesc

from AaronTools.atoms import BondOrder
from AaronTools.const import RADII

from scipy.spatial.distance import pdist

import numpy as np

from ora_stuff.mouse_modes import (
    SetHalfBond,
    SetSingleBond,
    SetAromaticBond,
    SetDoubleBond,
    SetTripleBond,
)

guessBondOrders_description = CmdDesc(
    required=[("selection", AtomsArg)],
    keyword=[
        ("coordinateSet", Or(IntArg, EnumOf(["all"]))),
        ("bondTolerance", FloatArg),
        ("TSBondTolerance", FloatArg),
    ],
    synopsis="guess bond orders based on the coordinates",
)

def guessBondOrders(
    session,
    selection,
    coordinateSet=None,
    bondTolerance=0.35,
    TSBondTolerance=0.6,
):
    """
    draw a TS bond
    selection - atoms
    coordinateSet - cd_id to use when guessing bond orders. If not given, current cs_id is used
    """

    bo_data = BondOrder()

    for structure, atoms in selection.by_structure:
        if coordinateSet is None:
            coordset_ids = [structure.active_coordset_id]
        elif coordinateSet == "all":
            coordset_ids = structure.coordset_ids
        else:
            try:
                coordset_ids = [coordinateSet]
                if int(coordinateSet) not in structure.coordset_ids:
                    session.logger.error("coordinateSet %s is not a valid coordinate set for %s" % (
                        coordinateSet,
                        structure.atomspec,
                    ))
                    return
            except ValueError:
                session.logger.error("expected 'all' or an integer for coordinateSet, got %s" % coordinateSet)
                return
        
        max_connected = dict()
        max_ts_distance = dict()
        for i, a1 in enumerate(atoms):
            for a2 in atoms[:i + 1]:
                key = bo_data.key(a1.element.name, a2.element.name)
                try:
                    max_connected[key]
                except KeyError:
                    max_connected[key] = (
                        RADII[a1.element.name] + RADII[a2.element.name] + bondTolerance
                    ) ** 2
                    max_ts_distance[key] = (
                        RADII[a1.element.name] + RADII[a2.element.name] + TSBondTolerance
                    ) ** 2
        
        selected = [structure.atoms.index(a) for a in atoms]
        
        for cs_id in coordset_ids:
            if coordinateSet is None:
                set_pb_coordsets = structure.coordset_ids
            elif coordinateSet == "all":
                set_pb_coordsets = [cs_id]
            else:
                set_pb_coordsets = structure.coordset_ids
        
            coords = structure.coordset(cs_id).xyzs
            distances = pdist(coords[selected], "sqeuclidean")
            
            bonds = {
                "broken": [],
                "half": [],
                "single": [],
                "aromatic": [],
                "double": [],
                "triple": [],
            }
            
            for i, a1 in enumerate(atoms):
                for j, a2 in enumerate(atoms[:i + 1]):
                    if i == j:
                        continue
                    key = bo_data.key(a1.element.name, a2.element.name)
                    # i, j = sorted([i, j])
                    ndx = len(atoms) * j + i - ((j + 2) * (j + 1)) // 2
                    # not bonded
                    if distances[ndx] > max_ts_distance[key] and TSBondTolerance > bondTolerance:
                        bonds["broken"].append((a1, a2))
                        continue
                    # ts bonded
                    if distances[ndx] < max_ts_distance[key] and distances[ndx] > max_connected[key]:
                        bonds["half"].append((a1, a2))
                        continue
                    # bonded and one of the atoms is H - can only have single bonds
                    if a1.element.name == "H" or a2.element.name == "H":
                        bonds["single"].append((a1, a2))
                        continue
                    
                    # bonded and order is guessed based on distance
                    # unless we don't have data for this pair of elements
                    try:
                        possible_bond_orders = bo_data.bonds[key]
                    except KeyError:
                        bonds["single"].append((a1, a2))
                        continue
                        
                    d = np.sqrt(distances[ndx])
                    closest = (0, None)
                    for order, length in possible_bond_orders.items():
                        diff = abs(length - d)
                        if closest[1] is None or diff < closest[1]:
                            closest = (order, diff)
                    
                    if closest[0] == "1.0":
                        bonds["single"].append((a1, a2))
                    elif closest[0] == "1.5":
                        bonds["aromatic"].append((a1, a2))
                    elif closest[0] == "2.0":
                        bonds["double"].append((a1, a2))
                    elif closest[0] == "3.0":
                        bonds["triple"].append((a1, a2))
                    else:
                        session.logger.warning("unanticipated bond order: %s. Expected one of %s" % (
                            closest[0], "1.0, 1.5, 2.0, 3.0",
                        ))
                    
            for bo in bonds:
                if not bonds[bo]:
                    continue
                
                if bo != "broken":
                    pbg = structure.pseudobond_group(
                        bo,
                        create_type=2,
                    )
    
                for (a1, a2) in bonds[bo]:
                    for bond in structure.bonds:
                        if a1 in bond.atoms and a2 in bond.atoms:
                            bond.delete()
                            break
                
                    if bo != "broken":
                        for cs in set_pb_coordsets:
                            if coordinateSet == "all" and cs < cs_id:
                                continue
                            pbg.new_pseudobond(a1, a2, cs_id=cs)

                            for other_bo in bonds:
                                if other_bo == bo:
                                    continue
                                other_pbg = structure.pseudobond_group(
                                    other_bo,
                                    create_type=None,
                                )
                                if not other_pbg:
                                    continue

                                pbs = other_pbg.get_pseudobonds(cs)
                                for pb in pbs:
                                    if a1 in pb.atoms and a2 in pb.atoms:
                                        pb.delete()
                
                if bo == "half":
                    pbg.color = SetHalfBond.color
                    pbg.dashes = SetHalfBond.dashes
                elif bo == "single":
                    pbg.color = SetSingleBond.color
                    pbg.dashes = SetSingleBond.dashes
                elif bo == "aromatic":
                    pbg.color = SetAromaticBond.color
                    pbg.dashes = SetAromaticBond.dashes
                elif bo == "double":
                    pbg.color = SetDoubleBond.color
                    pbg.dashes = SetDoubleBond.dashes
                elif bo == "triple":
                    pbg.color = SetTripleBond.color
                    pbg.dashes = SetTripleBond.dashes
