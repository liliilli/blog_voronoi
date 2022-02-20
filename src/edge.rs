use std::{cell::RefCell, rc::Rc};

use nalgebra::{Matrix3, Point2, Point3, Rotation2, Vector2};
type FPoint2 = Point2<f32>;
type FPoint3 = Point3<f32>;
type FVector2 = Vector2<f32>;
type FMatrix3 = Matrix3<f32>;

pub fn is_nearly_same_fpoint2(lhs: FPoint2, rhs: FPoint2, epsilon: f32) -> bool {
    let x = lhs[0] - rhs[0];
    let y = lhs[1] - rhs[1];

    (x.abs() <= epsilon) && (y.abs() <= epsilon)
}

#[derive(Clone, Copy, Debug)]
pub struct Edge {
    pub start: FPoint2,
    pub end: FPoint2,
}

impl Edge {
    pub fn new(start: FPoint2, end: FPoint2) -> Self {
        Edge { start, end }
    }

    pub fn get_offset(&self) -> FVector2 {
        self.end - self.start
    }

    pub fn is_nearly_same(&self, other: &Self) -> bool {
        // startとendが逆の可能性もあるので、それも踏まえてチェックしていく。
        // 最初はendまでのoffsetで1（並行）か-1(反転)でedgeが同じ位置に位置しているかを確認。
        let dir = self.get_offset().normalize();
        let other_dir = other.get_offset().normalize();
        let cosrad = dir.dot(&other_dir);

        if cosrad >= (1f32 - f32::EPSILON) {
            // 並行
            let is_start_nearly_equal =
                is_nearly_same_fpoint2(self.start, other.start, f32::EPSILON);
            let is_end_nearly_equal = is_nearly_same_fpoint2(self.end, other.end, f32::EPSILON);
            return is_start_nearly_equal && is_end_nearly_equal;
        } else if cosrad <= (-1f32 + f32::EPSILON) {
            // 反転
            let is_start_nearly_equal = is_nearly_same_fpoint2(self.start, other.end, f32::EPSILON);
            let is_end_nearly_equal = is_nearly_same_fpoint2(self.end, other.start, f32::EPSILON);
            return is_start_nearly_equal && is_end_nearly_equal;
        }

        false
    }

    pub fn get_reversed(&self) -> Self {
        Self {
            start: self.end,
            end: self.start,
        }
    }

    pub fn get_center(&self) -> FPoint2 {
        self.start + ((self.end - self.start) * 0.5f32)
    }

    pub fn try_get_direction(&self) -> Option<FVector2> {
        let offset = self.get_offset();
        if offset.magnitude_squared().is_normal() == false {
            return None;
        }

        Some(offset.normalize())
    }

    pub fn try_get_coefficients(&self) -> Option<(f32, f32, f32)> {
        // もしoffsetが0に近い場合には、coefficientが取れない。
        // またoffsetはyの符号を変える。
        try_get_coefficients(&(self.end - self.start), &self.start)
    }
}

pub fn try_get_coefficients(offset: &FVector2, point: &FPoint2) -> Option<(f32, f32, f32)> {
    // もしoffsetが0に近い場合には、coefficientが取れない。
    // またoffsetはyの符号を変える。
    let offset = offset.component_mul(&FVector2::new(1f32, -1f32));
    if offset.magnitude_squared() < f32::EPSILON {
        return None;
    }

    // offsetのどちらかで除算し、scaleをかけるが
    // xまたはyが0の可能性もあるし、なのでabsした値から大きいもので除算する。
    let unscaled_c = offset.dot(&(point - FPoint2::origin()));
    let factor = offset[0].abs().max(offset[1].abs());

    // scaleして返す。
    let scaled_ab = offset.scale(factor);
    Some((scaled_ab[0], scaled_ab[1], -1f32 * unscaled_c / factor))
}

#[derive(Debug)]
pub enum VoronoiEdge {
    None,
    InfinityFrom(FPoint2),
    Closed(Edge),
}

#[derive(Debug)]
pub struct EdgeContainer {
    pub site_edge: Edge, // start is bottom, end is top when not reversed.
    voronoi_edge: VoronoiEdge,
    pub bisector_pos: FPoint2,
    pub bisector_dir: FVector2,
}
pub type EdgeContainerRcCell = Rc<RefCell<EdgeContainer>>;

impl EdgeContainer {
    pub fn new(edge: &Edge) -> Self {
        // Check edge should be reversed or not.
        //
        // If start to edge offset does not satisfy below conditions, should be reversed.
        // * y-offset is 0, x-offset should be positive.
        // * x-offset is 0, y-offset should be positive.
        // * x-offset is not 0, y-offset should be positive.
        //
        // This logic is required for checking intersection, or left/right side for given point.
        // endはいつもstartよりxが大きいか、yが大きいことになる。
        let site_edge = {
            let to_end = edge.end - edge.start;
            let should_be_reversed = {
                // offset must not be Inf, or Nan.
                if to_end[1].is_normal() {
                    to_end[1] > 0f32
                } else {
                    to_end[0] > 0f32
                }
            };

            match should_be_reversed {
                true => *edge,
                false => edge.get_reversed(),
            }
        };

        // bisectorの方向を取得する。
        let bisector_pos = site_edge.get_center();
        let bisector_dir = {
            let d = site_edge.try_get_direction().unwrap();
            if site_edge.end[0] > site_edge.start[0] {
                // CCW
                Rotation2::new(90f32.to_radians()) * d
            } else {
                // CW
                Rotation2::new(-90f32.to_radians()) * d
            }
        };

        Self {
            site_edge,
            voronoi_edge: VoronoiEdge::None,
            bisector_pos,
            bisector_dir,
        }
    }

    pub fn into_rccell(self) -> EdgeContainerRcCell {
        Rc::new(RefCell::new(self))
    }

    pub fn is_left_of_bisect(&self, point: &FPoint2) -> Option<bool> {
        let bc = self.bisector_pos;
        let bd = self.bisector_dir;

        // If distance is negative, bisector is right side of point,
        // If distance is positive, bisector is left side of point.
        // If distance is 0, point is on the bisector.
        match try_get_coefficients(&bd, &bc) {
            Some((a, b, c)) => Some(((a * point[0]) + (b * point[1]) + c) < 0f32),
            None => None,
        }
    }
}
