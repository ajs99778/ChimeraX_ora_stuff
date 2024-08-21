from chimerax.atomic import AtomicStructureArg, AtomicStructuresArg
from chimerax.core.commands import IntArg, CmdDesc, run

import numpy as np

trim_cs_description = CmdDesc(
    required=[("selection", AtomicStructureArg)],
    keyword=[
        ("first", IntArg),
        ("last", IntArg),
    ],
    synopsis="delete a number of frames form the beginning or end of a trajectory",
)

flip_cs_description = CmdDesc(
    required=[("selection", AtomicStructureArg)],
    synopsis="reverse a trajectory",
)

combine_cs_description = CmdDesc(
    required=[
        ("selection", AtomicStructuresArg),
    ],
    synopsis="combine coordinate sets into a new structure",
)

def flip_cs(
    session,
    selection,
):
    current_coords = []
    for cs_id in selection.coordset_ids:
        current_coords.append(selection.coordset(cs_id).xyzs)
    current_coords.reverse()
    coords = np.array(current_coords)
    selection.add_coordsets(coords, replace=True)

def trim_cs(
    session,
    selection,
    first=0,
    last=0,
):
    current_coords = []
    for cs_id in selection.coordset_ids[first : selection.num_coordsets - last]:
        current_coords.append(selection.coordset(cs_id).xyzs)
    coords = np.array(current_coords)
    selection.add_coordsets(coords, replace=True)

def combine_cs(
    session,
    selection,
):
    if not all([selection[0].num_atoms == m.num_atoms for m in selection[1:]]):
        session.logger.error("cannot combine models: different numbers of atoms")
        return
    for i in range(0, selection[0].num_atoms):
        if not all([selection[0].atoms[i].element.name == m.atoms[i].element.name for m in selection[1:]]):
            session.logger.error("cannot combine models: atoms in a different order or molecules are different")
            return
    current_coords = []
    for m in selection:
        for cs_id in m.coordset_ids:
            current_coords.append(m.coordset(cs_id).xyzs)
    coords = np.array(current_coords)
    m = run(session, "combine %s" % selection[0].atomspec)
    run(session, "mcopy %s toAtoms %s" % (selection[0].atomspec, m.atomspec))
    m.add_coordsets(coords, replace=True)
