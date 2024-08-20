from chimerax.mouse_modes import MouseMode
from chimerax.atomic.structure import (
    PickedAtoms,
    PickedBonds,
    PickedAtom,
    PickedBond,
    PickedPseudobond,
)

from SEQCROW.mouse_modes import DrawBondMouseMode

bo_names = [
    'broken',
    'half',
    'single',
    'aromatic',
    'double',
    'triple',
]

class SetBondOrder(DrawBondMouseMode):
    name = None
    color = None
    dashes = 0
    draw_new = True

    def vr_press(self, event):
        self.mouse_up(event)

    def mouse_up(self, event):
        atom = False
        MouseMode.mouse_up(self, event)
        if self._atom1 is not None and self._atom1.deleted:
            self._atom1 = None

        x, y = event.position()
        pick = self.view.picked_object(x, y)
        
        if isinstance(pick, PickedPseudobond):
            bond = pick.pbond
            if bond.group.name == self.name:
                # run(self.session, "bond %s %s" % (bond.atoms[0].atomspec, bond.atoms[1].atomspec))
                bond.delete()
                self.reset()
            else:
                self.draw_new_pbond(*bond.atoms)
            return

        if isinstance(pick, PickedBond):
            bond = pick.bond
            self.draw_new_pbond(*bond.atoms)
            bond.delete()
            self.reset()
            return

        if isinstance(pick, PickedAtom) and not self._atom1:
            self._atom1 = pick.atom
            self.session.logger.status("drawing %s bond between %s and..." % (self.name, self._atom1.atomspec))
            return

        if (
            isinstance(pick, PickedAtom) and
            pick.atom.structure is self._atom1.structure and
            pick.atom is not self._atom1
        ):
            atom = pick.atom
        
        if not atom:
            self.session.logger.status("not atom clicked - not drawing a new %s bond" % self.name)
            return

        if atom.structure is self._atom1.structure and atom is not self._atom1:
            self.draw_new_pbond(self._atom1, atom)
            self.session.logger.status("drew new bond between %s and %s" % (
                atom.atomspec, self._atom1.atomspec,
            ))
            for bond in self._atom1.bonds:
                if atom in bond.atoms:
                    bond.delete()
                    break
            self.reset()

    def draw_new_pbond(self, atom1, atom2):
        if self.draw_new:
            pbg = atom1.structure.pseudobond_group(
                self.name,
                create_type=2
            )
            
            for cs_id in range(atom1.structure.active_coordset_id, atom1.structure.coordset_ids[-1] + 1):
                pbg.new_pseudobond(atom1, atom2, cs_id=cs_id)
            # bug in older versions of ChimeraX where new pseudobonds aren't 
            # displayed until something else changes
            # change the number of dashes
            pbg.dashes = self.dashes
            pbg.dashes += 2
            pbg.dashes -= 2
            pbg.color = self.color
        
        for kind in bo_names:
            if kind == self.name:
                continue
            other_pbg = atom1.structure.pseudobond_group(kind, create_type=None)
            if not other_pbg:
                continue
            for cs_id in range(atom1.structure.active_coordset_id, atom1.structure.coordset_ids[-1] + 1):
                pbs = other_pbg.get_pseudobonds(cs_id)
                for pb in pbs:
                    if atom1 in pb.atoms and atom2 in pb.atoms:
                        pb.delete()
            

class SetBrokenBond(SetBondOrder):
    name = "broken"
    color = [0, 0, 0, 1]
    dashes = 0
    draw_new = False


class SetHalfBond(SetBondOrder):
    name = "half"
    color = [50, 190, 50, 100]
    dashes = 0


class SetSingleBond(SetBondOrder):
    name = "single"
    color = [255, 255, 255, 255]
    dashes = 0


class SetAromaticBond(SetBondOrder):
    name = "aromatic"
    color = [255, 0, 255, 255]
    dashes = 4


class SetDoubleBond(SetBondOrder):
    name = "double"
    color = [0, 0, 0, 255]
    dashes = 4


class SetTripleBond(SetBondOrder):
    name = "triple"
    dashes = 6
    color = [255, 255, 0, 255]

