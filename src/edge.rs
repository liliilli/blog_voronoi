use std::{cell::RefCell, rc::Rc};

use itertools::Itertools;
use nalgebra::{ComplexField, Point2, Point3, Rotation2, Vector2};
type FPoint2 = Point2<f32>;
type FPoint3 = Point3<f32>;
type FVector2 = Vector2<f32>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Left,
    Up,
    Right,
    Down,
}

#[derive(Debug, Clone, Copy)]
pub enum VoronoiEdgeType {
    Closed,
    Opened(((Direction, FPoint2), bool /* isReversed */)),
    Boundary,
}

pub fn is_nearly_same_fpoint2(lhs: FPoint2, rhs: FPoint2, epsilon: f32) -> bool {
    let x = lhs[0] - rhs[0];
    let y = lhs[1] - rhs[1];

    (x.abs() <= epsilon) && (y.abs() <= epsilon)
}

static mut SITE_ID_COUNTER: usize = 0;

#[derive(Debug)]
pub struct Site {
    id: usize,
    pub point: FPoint2,
    pub voronoi_edges: Vec<(Edge, VoronoiEdgeType)>,
}
pub type SiteRcCell = Rc<RefCell<Site>>;

impl Site {
    pub fn into_rccell(self) -> SiteRcCell {
        Rc::new(RefCell::new(self))
    }

    pub fn from_point(point: FPoint2) -> Self {
        Self {
            id: unsafe {
                let id = SITE_ID_COUNTER;
                SITE_ID_COUNTER += 1;
                id
            },
            point,
            voronoi_edges: vec![],
        }
    }

    pub fn insert_edge(&mut self, mut edge: Edge, ve_type: VoronoiEdgeType) {
        if edge.start[0].abs() < f32::EPSILON {
            edge.start[0] = 0f32;
        }
        if edge.start[1].abs() < f32::EPSILON {
            edge.start[1] = 0f32;
        }
        if edge.end[0].abs() < f32::EPSILON {
            edge.end[0] = 0f32;
        }
        if edge.end[1].abs() < f32::EPSILON {
            edge.end[1] = 0f32;
        }
        self.voronoi_edges.push((edge, ve_type));
    }

    /// Check cell of site is closed or opened. 
    /// If site does not have any voronoi edges, just return [None] value.
    pub fn is_closed(&self) -> Option<bool> {
        if self.voronoi_edges.is_empty() {
            return None;
        }

        Some(self.voronoi_edges.iter().all(|(_, ve_type)| match ve_type {
            VoronoiEdgeType::Closed => true,
            // Boundary????????????
            _ => false,
        }))
    }

    pub fn find_open_edges_of_direction(&self, dir: Direction) -> Vec<(Edge, VoronoiEdgeType)>
    {
        if self.voronoi_edges.is_empty() {
            return vec![];
        }

        self.voronoi_edges.iter().filter(|(_, ve)| {
            match ve {
                VoronoiEdgeType::Opened(((target_dir, _), _)) => *target_dir == dir,
                _ => false,
            }
        }).map(|&v| v).collect_vec()
    }
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
        // start???end????????????????????????????????????????????????????????????????????????????????????
        // ?????????end?????????offset???1???????????????-1(??????)???edge???????????????????????????????????????????????????
        let dir = self.get_offset().normalize();
        let other_dir = other.get_offset().normalize();
        let cosrad = dir.dot(&other_dir);

        if cosrad >= (1f32 - f32::EPSILON) {
            // ??????
            let is_start_nearly_equal =
                is_nearly_same_fpoint2(self.start, other.start, f32::EPSILON);
            let is_end_nearly_equal = is_nearly_same_fpoint2(self.end, other.end, f32::EPSILON);
            return is_start_nearly_equal && is_end_nearly_equal;
        } else if cosrad <= (-1f32 + f32::EPSILON) {
            // ??????
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
        // ??????offset???0????????????????????????coefficient??????????????????
        // ??????offset???y????????????????????????
        try_get_coefficients(&(self.end - self.start), &self.start)
    }

}

pub fn try_get_coefficients(offset: &FVector2, point: &FPoint2) -> Option<(f32, f32, f32)> {
    // ??????offset???0????????????????????????coefficient??????????????????
    // ??????offset???y????????????????????????
    //let offset = //offset.component_mul(&FVector2::new(1f32, 1f32));
    if offset.magnitude_squared() < f32::EPSILON {
        return None;
    }

    // offset??????????????????????????????scale???????????????
    // x?????????y???0????????????????????????????????????abs????????????????????????????????????????????????
    let ab = offset.yx().component_mul(&FVector2::new(1f32, -1f32));
    let unscaled_c = ab.dot(&(point - FPoint2::origin()));
    let factor = ab[0].abs().max(ab[1].abs());

    // scale???????????????
    let scaled_ab = ab.scale(factor.recip());
    let mut a = scaled_ab[0];
    let mut b = scaled_ab[1];
    let mut c = -1f32 * unscaled_c / factor;
    if a.abs() < f32::EPSILON {
        a = 0f32;
    }
    if b.abs() < f32::EPSILON {
        b = 0f32;
    }
    if c.abs() < f32::EPSILON {
        c = 0f32;
    }

    // If b is negative, revert signnesses.
    if b.is_sign_negative() {
        a *= -1f32;
        b *= -1f32;
        c *= -1f32;
    }

    Some((a, b, c))
}

#[derive(Debug, Clone)]
pub struct SiteEdge {
    start: SiteRcCell,
    end: SiteRcCell,
}

impl SiteEdge {
    pub fn new(start: SiteRcCell, end: SiteRcCell) -> Self {
        Self { start, end }
    }

    pub fn start_point(&self) -> FPoint2 {
        self.start.borrow().point
    }

    pub fn end_point(&self) -> FPoint2 {
        self.end.borrow().point
    }

    pub fn get_center(&self) -> FPoint2 {
        let start = self.start_point();
        let end = self.end_point();

        start + ((end - start) * 0.5f32)
    }

    pub fn get_reversed(&self) -> Self {
        Self {
            start: self.end.clone(),
            end: self.start.clone(),
        }
    }

    pub fn get_offset(&self) -> FVector2 {
        let start = self.start_point();
        let end = self.end_point();

        end - start
    }

    pub fn try_get_direction(&self) -> Option<FVector2> {
        let offset = self.get_offset();
        if offset.magnitude_squared().is_normal() == false {
            return None;
        }

        Some(offset.normalize())
    }

    pub fn clone_start(&self) -> SiteRcCell {
        self.start.clone()
    }

    pub fn clone_end(&self) -> SiteRcCell {
        self.end.clone()
    }

    pub fn insert_edge(&mut self, edge: Edge, ve_type: VoronoiEdgeType) {
        self.start.borrow_mut().insert_edge(edge, ve_type);

        // ??????Opened?????????reversed??????????????????start???????????????????????????????????????
        let ve_type = match ve_type {
            VoronoiEdgeType::Opened((dir, rev)) => VoronoiEdgeType::Opened((dir, !rev)),
            v => v,
        };
        self.end.borrow_mut().insert_edge(edge, ve_type);
    }
}

#[derive(Debug, Clone)]
pub enum VoronoiEdge {
    None,
    InfFrom(Option<FPoint2>, Option<FPoint2>),
    Closed(Edge),
}

#[derive(Debug)]
pub struct EdgeContainer {
    pub site_edge: SiteEdge, // start is bottom, end is top when not reversed.
    pub voronoi_edge: VoronoiEdge,
    pub bisector_pos: FPoint2,
    pub bisector_dir: FVector2,
}
pub type EdgeContainerRcCell = Rc<RefCell<EdgeContainer>>;

impl EdgeContainer {
    pub fn new(in_edge: SiteEdge) -> Self {
        // Check edge should be reversed or not.
        //
        // If start to edge offset does not satisfy below conditions, should be reversed.
        // * y-offset is 0, x-offset should be positive.
        // * x-offset is 0, y-offset should be positive.
        // * x-offset is not 0, y-offset should be positive.
        //
        // This logic is required for checking intersection, or left/right side for given point.
        // end????????????start??????x??????????????????y??????????????????????????????
        let site_edge = {
            let to_end = in_edge.end_point() - in_edge.start_point();
            let should_be_reversed = {
                // offset must not be Inf, or Nan.
                if to_end[1].is_normal() {
                    to_end[1] > 0f32
                } else {
                    to_end[0] > 0f32
                }
            };

            match should_be_reversed {
                true => in_edge,
                false => in_edge.get_reversed(),
            }
        };

        // bisector???????????????????????????
        let bisector_pos = site_edge.get_center();
        let bisector_dir = {
            let d = site_edge.try_get_direction().unwrap();
            let mut dir = if site_edge.end_point()[0] > site_edge.start_point()[0] {
                // CCW
                Rotation2::new(90f32.to_radians()) * d
            } else {
                // CW
                Rotation2::new(-90f32.to_radians()) * d
            };

            // Check direction element is so small so errornous.
            if dir[0].abs() < f32::EPSILON {
                dir[0] = 0f32;
            }
            if dir[1].abs() < f32::EPSILON {
                dir[1] = 0f32;
            }
            dir
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

    pub fn is_bisect_left_of(&self, point: &FPoint2) -> Option<bool> {
        let bc = self.bisector_pos;
        let bd = self.bisector_dir;

        // If distance is negative, bisector is right side of point,
        // If distance is positive, bisector is left side of point.
        // If distance is 0, point is on the bisector.
        match try_get_coefficients(&bd, &bc) {
            Some((a, b, c)) => Some(((a * point[0]) + (b * point[1]) - c) > 0f32),
            None => None,
        }
    }

    pub fn try_get_coefficients_of_bisect(&self) -> Option<(f32, f32, f32)> {
        let bc = self.bisector_pos;
        let bd = self.bisector_dir;
        //print!("bc: {:?}, bd: {:?}", bc, bd);

        try_get_coefficients(&bd, &bc)
    }
}
