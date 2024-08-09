import os

from chimerax.atomic import AtomicStructureArg
from chimerax.core.toolshed import BundleAPI
from chimerax.core.commands import BoolArg, FileNameArg, FloatArg, StringArg

class _ora_stuff_API(BundleAPI):

    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        """
        custom initialization adds mouse modes to the menu
        """

        if session.ui.is_gui:
            from .mouse_modes import (
                SetBrokenBond,
                SetHalfBond,
                SetSingleBond,
                SetAromaticBond,
                SetDoubleBond,
                SetTripleBond,
            )
            
            session.ui.mouse_modes.add_mode(SetBrokenBond(session))
            session.ui.mouse_modes.add_mode(SetHalfBond(session))
            session.ui.mouse_modes.add_mode(SetSingleBond(session))
            session.ui.mouse_modes.add_mode(SetAromaticBond(session))
            session.ui.mouse_modes.add_mode(SetDoubleBond(session))
            session.ui.mouse_modes.add_mode(SetTripleBond(session))

    @staticmethod
    def run_provider(session, name, mgr, **kw):
        if mgr is session.presets:
            from .presets import run_preset
            run_preset(session, name, mgr)
        
        elif mgr is session.save_command:
            from chimerax.save_command import SaverInfo
            from .io import save_sdf, save_fbx
            
            if name == "SDF file":
                class Info(SaverInfo):
                    def save(self, session, path, **kw):
                        save_sdf(session, path, **kw)
                    
                    @property
                    def save_args(self):
                        return {
                            'model': AtomicStructureArg,
                            'tolerance': FloatArg,
                            'TStolerance': FloatArg,
                            'style': StringArg,
                            'coordsets': BoolArg,
                        }
                    
                    def save_args_widget(self, session):
                        from Qt.QtWidgets import (
                            QWidget, QFormLayout, QDoubleSpinBox, QCheckBox, QComboBox
                        )
                        from SEQCROW.widgets import ModelComboBox

                        widget = QWidget()
                        layout = QFormLayout(widget)
                        models = ModelComboBox(
                            session,
                            autoUpdate=False,
                        )
                        layout.addRow(models)
                        
                        tolerance = QDoubleSpinBox()
                        tolerance.setDecimals(2)
                        tolerance.setRange(0.0, 1.0)
                        tolerance.setSingleStep(0.05)
                        tolerance.setValue(0.35)
                        layout.addRow("bonded tolerance", tolerance)
                        
                        ts_tolerance = QDoubleSpinBox()
                        ts_tolerance.setDecimals(2)
                        ts_tolerance.setRange(0.0, 1.5)
                        ts_tolerance.setSingleStep(0.05)
                        ts_tolerance.setValue(0.6)
                        layout.addRow("half bond tolerance", ts_tolerance)
                        
                        style = QComboBox()
                        style.addItems([
                            "V2000",
                            "V3000",
                        ])
                        style.setCurrentIndex(1)
                        layout.addRow("style", style)
                        
                        coordsets = QCheckBox()
                        layout.addRow("all coordinate sets", coordsets)

                        return widget

                    def save_args_string_from_widget(self, widget):
                        from Qt.QtWidgets import QFormLayout
                        models = widget.layout().itemAt(0).widget().options_string()
                        tolerance = widget.layout().itemAt(1, QFormLayout.FieldRole).widget().value()
                        ts_tolerance = widget.layout().itemAt(2, QFormLayout.FieldRole).widget().value()
                        style = widget.layout().itemAt(3, QFormLayout.FieldRole).widget().currentText()
                        coordsets = widget.layout().itemAt(4, QFormLayout.FieldRole).widget().isChecked()
                        args = [
                            models.replace("models", "model", 1),
                            "tolerance", "%.2f" % tolerance,
                            "TStolerance", "%.2f" % ts_tolerance,
                            "style", style,
                            "coordsets", str(coordsets)
                        ]
                        return " ".join(args)

                return Info()

            if name == "FBX file":
                class Info(SaverInfo):
                    def save(self, session, path, **kw):
                        save_fbx(session, path, **kw)
                    
                    @property
                    def save_args(self):
                        return {
                            'model': AtomicStructureArg,
                            'blenderPath': FileNameArg,
                            'scriptOnly': BoolArg
                        }
                    
                    def save_args_widget(self, session):
                        from Qt.QtWidgets import QWidget, QFormLayout, QLineEdit, QCheckBox
                        from SEQCROW.widgets import ModelComboBox

                        widget = QWidget()
                        layout = QFormLayout(widget)
                        models = ModelComboBox(
                            session,
                            autoUpdate=False,
                        )
                        layout.addRow(models)
                        
                        # TODO: add a browse button
                        # I'm not sure if you can open a file select
                        # dialog while the file save dialog is open
                        blender = QLineEdit()
                        layout.addRow("blender executable", blender)
    
                        script_only = QCheckBox()
                        layout.addRow("python script only", script_only)
                        
                        return widget
    
                    def save_args_string_from_widget(self, widget):
                        from Qt.QtWidgets import QFormLayout
                        models = widget.layout().itemAt(0).widget().options_string()
                        blender = widget.layout().itemAt(1, QFormLayout.FieldRole).widget().text()
                        script_only = widget.layout().itemAt(2, QFormLayout.FieldRole).widget().isChecked()
                        args = [
                            models.replace("models", "model", 1),
                        ]
                        if blender and not script_only:
                            args.extend([
                                "blenderPath", '"' + blender + '"',
                                "scriptOnly", "false",
                            ])
                        else:
                            args.extend(["scriptOnly", "true"])
                        return " ".join(args)

                return Info()



bundle_api = _ora_stuff_API()
