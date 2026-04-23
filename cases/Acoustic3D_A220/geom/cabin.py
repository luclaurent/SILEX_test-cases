import gmsh
import numpy as np
from loguru import logger
from dataclasses import dataclass
import os
import tempfile
from typing import Any, cast

from seats import Seat, SeatPose, SeatShape

_LENGTH_CABIN = 15.0
_MINOR_AXIS = 3.28/2.0
_MAJOR_AXIS = 3.5/2.0
_Z_OFFSET_FLOOR = 2.13-3.5
_NB_LAYERS_EXTRUDE = 20


@dataclass
class CabinParameters:
    length: float = _LENGTH_CABIN
    half_width: float = _MINOR_AXIS
    half_height: float = _MAJOR_AXIS
    floor_z: float = _Z_OFFSET_FLOOR


@dataclass
class SeatLayoutParameters:
    seats_per_side: int = 3
    n_rows: int = 22
    row_pitch: float = 0.75
    front_margin: float = 1.7
    aisle_width: float = 0.51
    seat_base_z: float = _Z_OFFSET_FLOOR + 0.43
    side_clearance: float = 0.14
    include_overhead_bins: bool = True
    overhead_bin_height: float = 0.35
    overhead_bin_depth: float = 0.42
    overhead_bin_z: float = 0.88

class CabinInterior:
    """Parametric A220-like cabin interior model in Gmsh OCC."""

    def __init__(
        self,
        cabin_params: CabinParameters | None = None,
        seat_layout: SeatLayoutParameters | None = None,
        seat_shape: SeatShape | None = None,
    ):
        self.cabin = cabin_params or CabinParameters()
        self.seat_layout = seat_layout or SeatLayoutParameters()
        self.seat_shape = seat_shape or SeatShape()

        self.cabin_shell = []
        self.cabin_surfaces = []
        self.cabin_volume = None
        self.inside_air_volumes = []
        self.floor = None
        self.seats = []
        self.overhead_bins = []

    @property
    def length_cabin(self):
        return self.cabin.length

    @property
    def minor_axis(self):
        return self.cabin.half_width

    @property
    def major_axis(self):
        return self.cabin.half_height

    def create_cabin(self):
        if not gmsh.isInitialized():
            gmsh.initialize()

        x0 = -self.length_cabin / 2.0
        ellipse = gmsh.model.occ.addEllipse(x0, 0.0, 0.0, self.major_axis, self.minor_axis)

        # Rotate the section so width is on y and height on z.
        gmsh.model.occ.rotate([(1, ellipse)], x0, 0.0, 0.0, 0.0, 1.0, 0.0, np.pi / 2.0)

        ellipse_wire = gmsh.model.occ.addWire([ellipse])
        section = gmsh.model.occ.addPlaneSurface([ellipse_wire])
        extrusion = gmsh.model.occ.extrude(
            [(2, section)],
            self.length_cabin,
            0.0,
            0.0,
            numElements=[_NB_LAYERS_EXTRUDE],
            recombine=True,
        )

        self.cabin_shell = extrusion
        self.cabin_surfaces = [s[1] for s in extrusion if s[0] == 2]
        volume_list = [s[1] for s in extrusion if s[0] == 3]
        self.cabin_volume = volume_list[0] if volume_list else None
        gmsh.model.occ.synchronize()

    def create_floor(self):
        if self.cabin_volume is None:
            self.create_cabin()

        z_floor = self.cabin.floor_z
        floor_half_width = self._section_half_width_at_z(z_floor)
        self.floor = gmsh.model.occ.addRectangle(
            -self.length_cabin / 2.0,
            -floor_half_width,
            z_floor,
            self.length_cabin,
            2.0 * floor_half_width,
        )
        gmsh.model.occ.synchronize()

    def _section_half_width_at_z(self, z_value):
        z_ratio = z_value / self.major_axis
        if abs(z_ratio) >= 1.0:
            return 0.0
        return self.minor_axis * np.sqrt(1.0 - z_ratio**2)

    def create_overhead_bins(self):
        if not self.seat_layout.include_overhead_bins:
            return

        y_center = 0.5 * self.seat_layout.aisle_width + 1.5 * self.seat_shape.width
        y_center += self.seat_layout.side_clearance + 0.5 * self.seat_layout.overhead_bin_depth
        z0 = self.seat_layout.overhead_bin_z

        left = gmsh.model.occ.addBox(
            -self.length_cabin / 2.0,
            y_center - self.seat_layout.overhead_bin_depth,
            z0,
            self.length_cabin,
            self.seat_layout.overhead_bin_depth,
            self.seat_layout.overhead_bin_height,
        )
        right = gmsh.model.occ.addBox(
            -self.length_cabin / 2.0,
            -y_center,
            z0,
            self.length_cabin,
            self.seat_layout.overhead_bin_depth,
            self.seat_layout.overhead_bin_height,
        )
        self.overhead_bins = [(3, left), (3, right)]
        gmsh.model.occ.synchronize()

    def create_seats(self):
        seat_builder = Seat(gmsh, self.seat_shape)
        y_right = 0.5 * self.seat_layout.aisle_width + 0.5 * self.seat_shape.width
        y_left = -y_right

        for i_row in range(self.seat_layout.n_rows):
            x = -self.length_cabin / 2.0 + self.seat_layout.front_margin + i_row * self.seat_layout.row_pitch
            if x > self.length_cabin / 2.0 - 0.6:
                break

            for i_side in range(self.seat_layout.seats_per_side):
                y_pos_right = y_right + i_side * self.seat_shape.width
                y_pos_left = y_left - i_side * self.seat_shape.width

                self.seats.extend(
                    seat_builder.create(
                        SeatPose(
                            x=x,
                            y=y_pos_right,
                            z=self.seat_layout.seat_base_z,
                        )
                    )
                )
                self.seats.extend(
                    seat_builder.create(
                        SeatPose(
                            x=x,
                            y=y_pos_left,
                            z=self.seat_layout.seat_base_z,
                        )
                    )
                )

        gmsh.model.occ.synchronize()

    def build(self):
        self.create_cabin()
        self.create_floor()
        self.create_overhead_bins()
        self.create_seats()

        gmsh.model.addPhysicalGroup(3, [self.cabin_volume], name="CabinAirVolume")
        gmsh.model.addPhysicalGroup(2, [self.floor], name="CabinFloor")

        seat_volume_ids = [v[1] for v in self.seats if v[0] == 3]
        if seat_volume_ids:
            gmsh.model.addPhysicalGroup(3, seat_volume_ids, name="Seats")

        bin_volume_ids = [v[1] for v in self.overhead_bins if v[0] == 3]
        if bin_volume_ids:
            gmsh.model.addPhysicalGroup(3, bin_volume_ids, name="OverheadBins")

        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        return {
            "cabin_volume": self.cabin_volume,
            "floor_surface": self.floor,
            "seat_volumes": [v[1] for v in self.seats if v[0] == 3],
            "overhead_bin_volumes": [v[1] for v in self.overhead_bins if v[0] == 3],
        }

    def build_inside_air_domain(self, subtract_internal_solids=True):
        if self.cabin_volume is None:
            self.build()

        if self.inside_air_volumes:
            return self.inside_air_volumes

        solid_obstacles = [v for v in self.seats if v[0] == 3]
        solid_obstacles.extend(v for v in self.overhead_bins if v[0] == 3)

        if subtract_internal_solids and solid_obstacles:
            cut_result, _ = gmsh.model.occ.cut(
                [(3, self.cabin_volume)],
                solid_obstacles,
                removeObject=True,
                removeTool=True,
            )
            self.inside_air_volumes = [v[1] for v in cut_result if v[0] == 3]
        else:
            self.inside_air_volumes = [v for v in [self.cabin_volume] if v is not None]
            if solid_obstacles:
                gmsh.model.occ.remove(solid_obstacles, recursive=True)

        self.seats = []
        self.overhead_bins = []

        if self.floor is not None:
            gmsh.model.occ.remove([(2, self.floor)], recursive=False)
            self.floor = None

        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()
        return self.inside_air_volumes

    def show_all_data(self):
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()
        txt = []
        txt.extend(self.show_data(gmsh.model, name='GEO'))
        txt.extend(self.show_data(gmsh.model.occ, name='OCC'))
        return txt

    def show_data(self, model=None, name='OCC'):
        if model is None:
            model = gmsh.model.occ
        
        list_points = [x[1] for x in model.getEntities(0)]
        list_curves = [x[1] for x in model.getEntities(1)]
        list_surfaces = [x[1] for x in model.getEntities(2)]
        list_volumes = [x[1] for x in model.getEntities(3)]
        txt = []
        if len(list_points)>0:
            txt.append(f"Points ({name}): {list_points}")
            logger.info(txt[-1])
        if len(list_curves)>0:
            txt.append(f"Curves ({name}): {list_curves}")
            logger.info(txt[-1])
        if len(list_surfaces)>0:
            txt.append(f"Surfaces ({name}): {list_surfaces}")
            logger.info(txt[-1])
        if len(list_volumes)>0:
            txt.append(f"Volumes ({name}): {list_volumes}")
            logger.info(txt[-1])
        return txt

    def export_mesh(self, filename="cabin.msh", mesh_size=0.22):
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.5 * mesh_size)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)
        gmsh.model.mesh.generate(3)
        gmsh.write(filename)
        logger.info(f"Wrote mesh to {filename}")

    def _prepare_mesh(self, mesh_size=0.22, inside_air=False, subtract_internal_solids=False):
        if inside_air:
            self.build_inside_air_domain(
                subtract_internal_solids=subtract_internal_solids
            )
        elif self.cabin_volume is None:
            self.build()

        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.5 * mesh_size)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)
        gmsh.model.mesh.generate(3)

    def show_geo_gmsh(self):
        if self.cabin_volume is None:
            self.build()
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()
        gmsh.fltk.run()

    def show_mesh_gmsh(self, mesh_size=0.22, inside_air=False, subtract_internal_solids=False):
        self._prepare_mesh(
            mesh_size=mesh_size,
            inside_air=inside_air,
            subtract_internal_solids=subtract_internal_solids,
        )
        gmsh.fltk.run()

    def _mesh_to_vtk_file(self, mesh_size=0.22, inside_air=False, subtract_internal_solids=False):
        self._prepare_mesh(
            mesh_size=mesh_size,
            inside_air=inside_air,
            subtract_internal_solids=subtract_internal_solids,
        )

        with tempfile.NamedTemporaryFile(suffix=".vtk", delete=False) as handle:
            vtk_path = handle.name
        gmsh.write(vtk_path)
        return vtk_path

    def show_mesh_pyvista(
        self,
        mesh_size=0.22,
        inside_air=False,
        subtract_internal_solids=False,
        show_edges=True,
        opacity=1.0,
    ):
        try:
            import pyvista as pv
        except ModuleNotFoundError as exc:
            raise RuntimeError("PyVista is required for show_mesh_pyvista().") from exc

        vtk_path = self._mesh_to_vtk_file(
            mesh_size=mesh_size,
            inside_air=inside_air,
            subtract_internal_solids=subtract_internal_solids,
        )
        try:
            grid = cast(Any, pv.read(vtk_path))
            pv.plot(grid, show_edges=show_edges, opacity=opacity)
        finally:
            if os.path.exists(vtk_path):
                os.remove(vtk_path)

    def show_geo_pyvista(self, surface_mesh_size=0.30, show_edges=True):
        try:
            import pyvista as pv
        except ModuleNotFoundError as exc:
            raise RuntimeError("PyVista is required for show_geo_pyvista().") from exc

        if self.cabin_volume is None:
            self.build()

        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.5 * surface_mesh_size)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", surface_mesh_size)
        gmsh.model.mesh.generate(2)

        with tempfile.NamedTemporaryFile(suffix=".vtk", delete=False) as handle:
            vtk_path = handle.name
        gmsh.write(vtk_path)

        try:
            mesh = cast(Any, pv.read(vtk_path))
            pv.plot(mesh, style="wireframe", show_edges=show_edges, color="black")
        finally:
            if os.path.exists(vtk_path):
                os.remove(vtk_path)

    def export_inside_air_mesh(
        self,
        filename="cabin_inside_air.msh",
        mesh_size=0.22,
        subtract_internal_solids=False,
    ):
        air_volumes = self.build_inside_air_domain(
            subtract_internal_solids=subtract_internal_solids
        )

        # Export only inside-air elements by keeping a single physical volume group.
        existing_groups = gmsh.model.getPhysicalGroups()
        if existing_groups:
            gmsh.model.removePhysicalGroups(existing_groups)
        gmsh.model.addPhysicalGroup(3, air_volumes, name="InsideAir")
        gmsh.option.setNumber("Mesh.SaveAll", 0)

        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.5 * mesh_size)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)
        gmsh.model.mesh.generate(3)
        gmsh.write(filename)

        nb_elements = 0
        for volume_tag in air_volumes:
            _, element_tags, _ = gmsh.model.mesh.getElements(3, volume_tag)
            nb_elements += sum(len(tags) for tags in element_tags)

        logger.info(
            f"Wrote inside-air mesh to {filename} with {len(air_volumes)} volumes and {nb_elements} elements"
        )

        return {
            "filename": filename,
            "inside_air_volumes": air_volumes,
            "n_elements": nb_elements,
        }
        
        
    def show(self):
        self.show_all_data()
        self.show_geo_gmsh()
        
if __name__ == "__main__":
    gmsh.model.add("A220_cabin_interior")
    interior = CabinInterior()
    interior.build()
    interior.show()