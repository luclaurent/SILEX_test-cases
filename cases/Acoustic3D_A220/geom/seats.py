from dataclasses import dataclass


@dataclass
class SeatShape:
    """Parametric dimensions of a simplified aircraft seat."""

    width: float = 0.43
    depth: float = 0.48
    cushion_thickness: float = 0.12
    back_thickness: float = 0.08
    back_height: float = 0.72
    leg_radius: float = 0.015
    leg_height: float = 0.43
    headrest_height: float = 0.16


@dataclass
class SeatPose:
    """Seat placement in the cabin frame."""

    x: float
    y: float
    z: float


class Seat:
    """Build a single seat as a set of Gmsh OCC volumes."""

    def __init__(self, gmsh_module, shape: SeatShape | None = None):
        self.gmsh = gmsh_module
        self.shape = shape or SeatShape()

    def create(self, pose: SeatPose):
        occ = self.gmsh.model.occ
        s = self.shape

        seat_volumes = []

        # Cushion block.
        cushion = occ.addBox(
            pose.x - 0.5 * s.depth,
            pose.y - 0.5 * s.width,
            pose.z,
            s.depth,
            s.width,
            s.cushion_thickness,
        )
        seat_volumes.append((3, cushion))

        # Back block.
        back = occ.addBox(
            pose.x + 0.5 * s.depth - s.back_thickness,
            pose.y - 0.5 * s.width,
            pose.z + s.cushion_thickness,
            s.back_thickness,
            s.width,
            s.back_height,
        )
        seat_volumes.append((3, back))

        # Headrest block.
        headrest = occ.addBox(
            pose.x + 0.5 * s.depth - s.back_thickness,
            pose.y - 0.5 * s.width,
            pose.z + s.cushion_thickness + s.back_height - s.headrest_height,
            s.back_thickness,
            s.width,
            s.headrest_height,
        )
        seat_volumes.append((3, headrest))

        # Two support legs.
        y_leg_offset = 0.3 * s.width
        x_leg = pose.x - 0.25 * s.depth
        for y_leg in (pose.y - y_leg_offset, pose.y + y_leg_offset):
            leg = occ.addCylinder(
                x_leg,
                y_leg,
                pose.z - s.leg_height,
                0.0,
                0.0,
                s.leg_height,
                s.leg_radius,
            )
            seat_volumes.append((3, leg))

        return seat_volumes
