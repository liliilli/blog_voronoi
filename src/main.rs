#![feature(drain_filter)]
use ordered_float::OrderedFloat;

use std::{
    cell::RefCell,
    collections::{BTreeMap, VecDeque, HashSet},
    rc::Rc,
};

use itertools::Itertools;
use nalgebra::{Point2, Vector2};
type FPoint2 = Point2<f32>;
type FVector2 = Vector2<f32>;

mod edge;
use edge::{Edge, EdgeContainer, EdgeContainerRcCell, VoronoiEdge};

///
type HalfEdgeRcCell = Rc<RefCell<HalfEdge>>;

/// KeyタイプはHalfEdgeの概略なx軸位置を示す。
/// 各Siteの一番近そうなHalfEdgeを取得するため。
type HalfEdgeBTreeMap = BTreeMap<OrderedFloat<f32>, HalfEdgeRcCell>;

///
type HalfEdgeVertexEventRcCell = Rc<RefCell<HalfEdgeVertexEvent>>;
type HalfEdgeVertexEventBTreeMap = BTreeMap<OrderedFloat<f32>, HalfEdgeVertexEventRcCell>;

static mut HE_ID_COUNTER: usize = 0;
static mut VE_ID_COUNTER: usize = 0;

#[derive(Default)]
struct HalfEdge {
    id: usize,
    left_halfedge: Option<HalfEdgeRcCell>,
    right_halfedge: Option<HalfEdgeRcCell>,
    pub edge: Option<EdgeContainerRcCell>,
    pub is_reversed_he: bool,
    pub ve_ref: Option<HalfEdgeVertexEventRcCell>,
}

impl std::fmt::Debug for HalfEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Just show left or right is exist.
        let left_id = match self.left_halfedge.as_ref() {
            Some(v) => Some(v.borrow().id),
            None => None,
        };
        let right_id = match self.right_halfedge.as_ref() {
            Some(v) => Some(v.borrow().id),
            None => None,
        };
        let ve_id = match self.ve_ref.as_ref() {
            Some(ve) => Some(ve.borrow().id),
            None => None,
        };
        let is_boundary = self.edge.is_none();

        f.debug_struct("HE")
            .field("id", &self.id)
            .field("boundary", &is_boundary)
            .field("edge", &self.edge)
            .field("is_reversed", &self.is_reversed_he)
            .field("lid", &left_id)
            .field("rid", &right_id)
            .field("veid", &ve_id)
            .finish()
    }
}

impl HalfEdge {
    pub fn default_rccell() -> HalfEdgeRcCell {
        Rc::new(RefCell::new(Self::new()))
    }

    pub fn into_rccell(self) -> HalfEdgeRcCell {
        Rc::new(RefCell::new(self))
    }

    fn new() -> Self {
        Self {
            id: unsafe {
                let id = HE_ID_COUNTER;
                HE_ID_COUNTER += 1;
                id
            },
            ..Default::default()
        }
    }

    pub fn from_edge(edge: &Edge, is_reversed_he: bool) -> Self {
        Self {
            id: unsafe {
                let id = HE_ID_COUNTER;
                HE_ID_COUNTER += 1;
                id
            },
            edge: Some(EdgeContainer::new(edge).into_rccell()),
            is_reversed_he,
            ..Default::default()
        }
    }

    pub fn from_edge_container(ec: EdgeContainerRcCell, is_reversed_he: bool) -> Self {
        Self {
            id: unsafe {
                let id = HE_ID_COUNTER;
                HE_ID_COUNTER += 1;
                id
            },
            edge: Some(ec),
            is_reversed_he,
            ..Default::default()
        }
    }

    pub fn clone_edge(&self) -> Option<EdgeContainerRcCell> {
        match self.edge.as_ref() {
            Some(ec) => Some(ec.clone()),
            None => None,
        }
    }

    pub fn chain_as_right(left: &mut HalfEdgeRcCell, right: HalfEdgeRcCell) {
        // Get old right from left
        let left_old_right = match &left.borrow_mut().right_halfedge {
            Some(or) => Some(or.clone()),
            None => None,
        };

        // Connect to each other. (Left->Right)
        left.borrow_mut().right_halfedge = Some(right.clone());

        {
            // もしRightのLeftがすでにあれば、LeftのRightの連結を解除する。
            let mut right_mut = right.borrow_mut();
            match &mut right_mut.left_halfedge {
                Some(ol) => ol.borrow_mut().right_halfedge = None,
                _ => (),
            };

            // Connect to each other. (Right->Left)
            right_mut.left_halfedge = Some(left.clone());
        }

        // もしOld-Rightが存在していれば、
        // OldRight->left = right, right->right = old-rightにする。
        if let Some(right_right) = left_old_right {
            right_right.borrow_mut().left_halfedge = Some(right.clone());
            right.borrow_mut().right_halfedge = Some(right_right);
        }
    }

    pub fn is_bisect_left_of(&self, point: &FPoint2) -> Option<bool> {
        match &self.edge {
            Some(container) => {
                // 単純にbisectから左なのかを計算するのではなく、
                // reversedかではないかによって特殊に判定を行う必要がある。
                let c = container.borrow();
                let topsite = &c.site_edge.end;

                let is_right_of_site = point[0] > topsite[0];
                if is_right_of_site && !self.is_reversed_he {
                    return Some(true);
                } else if !is_right_of_site && self.is_reversed_he {
                    return Some(false);
                }

                match c.is_bisect_left_of(point) {
                    Some(v) => Some(match self.is_reversed_he {
                        true => !v,
                        false => v,
                    }),
                    None => None,
                }
            }
            None => None,
        }
    }

    pub fn try_get_edge_start(&self) -> Option<FPoint2> {
        match self.edge.as_ref() {
            Some(ec) => {
                let borrowed_ec = ec.borrow();
                Some(match self.is_reversed_he {
                    false => borrowed_ec.site_edge.start,
                    true => borrowed_ec.site_edge.end,
                })
            }
            None => None,
        }
    }

    pub fn try_get_edge_end(&self) -> Option<FPoint2> {
        match self.edge.as_ref() {
            Some(ec) => {
                let borrowed_ec = ec.borrow();
                Some(match self.is_reversed_he {
                    true => borrowed_ec.site_edge.start,
                    false => borrowed_ec.site_edge.end,
                })
            }
            None => None,
        }
    }

    pub fn try_get_bisect_intersected(&self, rhs: &Self) -> Option<FPoint2> {
        let (lhs_e, lbv, lbd) = {
            match self.edge.as_ref() {
                Some(ec) => {
                    let borrowed = ec.borrow();
                    (
                        borrowed.site_edge,
                        borrowed.bisector_pos,
                        borrowed.bisector_dir,
                    )
                }
                None => return None,
            }
        };
        let (rhs_e, rbv, rbd) = {
            match rhs.edge.as_ref() {
                Some(ec) => {
                    let borrowed = ec.borrow();
                    (
                        borrowed.site_edge,
                        borrowed.bisector_pos,
                        borrowed.bisector_dir,
                    )
                }
                None => return None,
            }
        };

        // Check lhs_e and rhs_e has same parent.
        if lhs_e.end == rhs_e.end {
            return None;
        }

        // Check lhs_e and rhs_e is parallel to each other.
        // 元コードではperp dot productで0に近いかで確認してた。
        {
            let lhs_dir = lhs_e.try_get_direction().unwrap();
            let rhs_dir = rhs_e.try_get_direction().unwrap();
            if lhs_dir.dot(&rhs_dir).abs() >= (1f32 - f32::EPSILON) {
                return None;
            }
        }

        // その後に計算されるsとtは0..=lhs_lenと0..=rhs_lenの中に入らなければならない。
        // From real-time rendering, 4th ed.
        let s = (rbv - lbv).perp(&rbd) / lbd.perp(&rbd);
        let t = (lbv - rbv).perp(&lbd) / rbd.perp(&lbd);
        if !s.is_finite() || !t.is_finite() {
            return None;
        }
        let intersected = lbv + (lbd * s);

        // もしEdgeのend点が逆なら、ターゲットを逆にする。
        let (el, e) = {
            let le = &lhs_e.end;
            let re = &rhs_e.end;

            if le[1] < re[1] || (le[1] == re[1] && le[0] < re[0]) {
                (self, &lhs_e)
            } else {
                (rhs, &rhs_e)
            }
        };

        // Right of site of e?
        if intersected[0] >= e.end[0] {
            if !el.is_reversed_he {
                return None;
            }
        } else {
            if el.is_reversed_he {
                return None;
            }
        }

        Some(intersected)
    }

    pub fn update_voronoi_edge(&mut self, point: FPoint2, reversed: bool) -> Option<Edge> {
        // tt => f, tf = t, ff => f, ft => t.
        let is_reversed_he = self.is_reversed_he ^ reversed;
        match self.edge.as_mut() {
            Some(ec) => {
                let mut mut_ec = ec.borrow_mut();
                let next = match &mut_ec.voronoi_edge {
                    VoronoiEdge::None => match is_reversed_he {
                        true => VoronoiEdge::InfFrom(None, Some(point)),
                        false => VoronoiEdge::InfFrom(Some(point), None),
                    },
                    VoronoiEdge::InfFrom(s, None) => match is_reversed_he {
                        true => VoronoiEdge::Closed(Edge::new(s.unwrap(), point)),
                        false => unreachable!("Unreachable! InfFrom"),
                    },
                    &VoronoiEdge::InfFrom(None, e) => match is_reversed_he {
                        false => VoronoiEdge::Closed(Edge::new(point, e.unwrap())),
                        true => unreachable!("Unreachable! InfFrom"),
                    },
                    _ => unreachable!("Unreachable! Closed"),
                };
                mut_ec.voronoi_edge = next.clone();

                // もしedgeが完成されたら、返す。
                if let VoronoiEdge::Closed(e) = next {
                    return Some(e);
                }

                None
            }
            None => None,
        }
    }

    pub fn try_get_voronoi_edge(&self, min_boundary: FPoint2, max_boundary: FPoint2) -> Option<Edge> {
        let (mut p1, mut p2) = {
            let ec = self.edge.as_ref().unwrap();
            let borrowed_ec = ec.borrow();

            (borrowed_ec.site_edge.start, borrowed_ec.site_edge.end)
        };

        let (a, b, c) = {
            let ec = self.edge.as_ref().unwrap();
            let borrowed_ec = ec.borrow();
            borrowed_ec.try_get_coefficients_of_bisect().unwrap()
        };
        if !a.is_normal() && !b.is_normal() {
            return None;
        }

        let (s1, s2) = {
            let ec = self.edge.as_ref().unwrap();
            let borrowed_ec = ec.borrow();
            match borrowed_ec.voronoi_edge {
                VoronoiEdge::InfFrom(s, e) => {
                    (s.clone(), e.clone())
                },
                _ => unreachable!(),
            }
        };

        // Get a line.
        if b.is_normal() //a.is_subnormal()
        {
            p1 = {
                let mut x = min_boundary[0];
                if let Some(p) = s1 {
                    if p[0] > min_boundary[0]  {
                        x = p[0];
                    }
                }
                let x = x.clamp(min_boundary[0], max_boundary[0]);
                let y = ((a * x) + c) / b * -1f32;
                FPoint2::new(x, y)
            };
            p2 = {
                let mut x = max_boundary[0];
                if let Some(p) = s2 {
                    if p[0] < max_boundary[0] {
                        x = p[0];
                    }
                }
                let x = x.clamp(min_boundary[0], max_boundary[0]);
                let y = ((a * x) + c) / b * -1f32;
                FPoint2::new(x, y)
            };

            // Check updated (p1, p2) is out of bound.
            let (miny, maxy) = (min_boundary[1], max_boundary[1]);
            if (p1[1] > maxy && p2[1] > maxy) || (p1[1] < miny && p2[1] < miny)
            {
                return None;
            }

            // Clip by y1->x1, and y2->x2.
            // a.is_subnormal()の場合、上の条件式で弾かれるのでaは有効な値があると仮定して勧めて良いかも。
            let p1y_clipped = p1[1].clamp(miny, maxy);
            if p1y_clipped != p1[1] {
                let x = ((b * p1y_clipped) + c) / a * -1f32;
                p1 = FPoint2::new(x, p1y_clipped);
            }

            let p2y_clipped = p2[1].clamp(miny, maxy);
            if p2y_clipped != p2[1] {
                let x = ((b * p2y_clipped) + c) / a * -1f32;
                p2 = FPoint2::new(x, p2y_clipped);
            }
        }
        else if a.is_normal()
        {
            p1 = {
                let mut y = min_boundary[1];
                if let Some(p) = s1 {
                    if p[1] > min_boundary[1]  {
                        y = p[1];
                    }
                }
                let y = y.clamp(min_boundary[1], max_boundary[1]);
                let x = ((b * y) + c) / a * -1f32;
                FPoint2::new(x, y)
            };
            p2 = {
                let mut y = min_boundary[1];
                if let Some(p) = s2 {
                    if p[1] < max_boundary[1]  {
                        y = p[1];
                    }
                }
                let y = y.clamp(min_boundary[1], max_boundary[1]);
                let x = ((b * y) + c) / a * -1f32;
                FPoint2::new(x, y)
            };

            // Check updated (p1, p2) is out of bound.
            let (minx, maxx) = (min_boundary[0], max_boundary[0]);
            if (p1[0] > maxx && p2[0] > maxx) || (p1[0] < minx && p2[0] < minx)
            {
                return None;
            }

            // Clip by x1->y1, and x2->y2.
            // a.is_subnormal()の場合、上の条件式で弾かれるのでaは有効な値があると仮定して勧めて良いかも。
            let p1x_clipped = p1[0].clamp(minx, maxx);
            if p1x_clipped != p1[0] {
                let y = ((a * p1x_clipped) + c) / b * -1f32;
                p1 = FPoint2::new(p1x_clipped, y);
            }

            let p2x_clipped = p2[0].clamp(minx, maxx);
            if p2x_clipped != p2[0] {
                let y = ((a * p2x_clipped) + c) / b * -1f32;
                p2 = FPoint2::new(p2x_clipped, y);
            }
        }

        Some(Edge::new(p1, p2))
    }
}

#[derive(Debug)]
struct HalfEdgeMap {
    map: HalfEdgeBTreeMap,
    most_left: HalfEdgeRcCell,
    most_right: HalfEdgeRcCell,
}

impl HalfEdgeMap {
    pub fn new() -> Self {
        let mut map = HalfEdgeBTreeMap::new();

        let mut most_left = HalfEdge::default_rccell();
        let most_right = HalfEdge::default_rccell();
        HalfEdge::chain_as_right(&mut most_left, most_right.clone());

        const MIN_LIMIT: f32 = -5000f32;
        const MAX_LIMIT: f32 = 5000f32;

        map.insert(OrderedFloat(MIN_LIMIT), most_left.clone());
        map.insert(OrderedFloat(MAX_LIMIT), most_right.clone());

        Self {
            map,
            most_left,
            most_right,
        }
    }

    pub fn print_all(&self) {
        self.map.iter().for_each(|(k, v)| {
            println!("{:?}", (k, v));
            if v.borrow().right_halfedge.is_some() {
                let mut cursor = v.borrow().right_halfedge.as_ref().unwrap().clone();
                loop {
                    println!("Next: {:?}", cursor);
                    if cursor.borrow().right_halfedge.is_some() {
                        let next = cursor.borrow().right_halfedge.as_ref().unwrap().clone();
                        cursor = next;
                    } else {
                        break;
                    }
                }
            }
        });
    }

    pub fn check_validation(&self) {
        if self.map.is_empty() {
            return;
        }

        // Check from start to end using right_halfedge.
        let mut cursor = match self.most_left.borrow().right_halfedge.as_ref() {
            Some(right) => right.clone(),
            None => return,
        };
        if cursor.borrow().edge.is_none() {
            return;
        }

        let mut next_start_site = {
            let borrowed_cursor = cursor.borrow();
            let is_reversed_he = borrowed_cursor.is_reversed_he;
            let site_edge = borrowed_cursor.edge.as_ref().unwrap().borrow().site_edge;

            match is_reversed_he {
                true => site_edge.end,
                false => site_edge.start,
            }
        };
        let goal_site = next_start_site;

        loop 
        {
            if cursor.borrow().edge.is_none() { // Is most-right?
                break;
            }

            let (start, end) = {
                let borrowed_cursor = cursor.borrow();
                let is_reversed_he = borrowed_cursor.is_reversed_he;
                let site_edge = borrowed_cursor.edge.as_ref().unwrap().borrow().site_edge;

                match is_reversed_he {
                    false => {
                        (site_edge.start, site_edge.end)
                    },
                    true => {
                        (site_edge.end, site_edge.start)
                    },
                }
            };

            assert!(start == next_start_site);
            next_start_site = end;

            // Update cursor to right next.
            {
                let next = cursor.borrow().right_halfedge.as_ref().unwrap().clone();
                cursor = next;
            }
        }

        assert!(next_start_site == goal_site);
    }

    pub fn get_nearest_left_of<'a>(&'a self, site: &FPoint2) -> Option<HalfEdgeRcCell> {
        if self.map.is_empty() {
            return None;
        }

        // 左で近そうなものを取得する。
        let x = site[0];
        let mut maybe_left = match self.map.range(..OrderedFloat(x)).next() {
            Some((_, v)) => v.clone(),
            None => return None,
        };

        // 取得したものが本当に左に位置しているかを確認する。
        let p_most_left = self.most_left.as_ptr();
        let p_most_right = self.most_right.as_ptr();
        let is_most_left = maybe_left.as_ptr() == p_most_left;
        let is_left_of_point = maybe_left.as_ptr() != p_most_right
            && maybe_left.borrow().is_bisect_left_of(&site).unwrap_or(true);

        if is_most_left || is_left_of_point {
            loop {
                let right = maybe_left.borrow().right_halfedge.as_ref().unwrap().clone();
                maybe_left = right;

                // rightが本当にsiteの右なら、ループを止めてrightのleftを一番近いleftとしてみなす。
                if maybe_left.as_ptr() == p_most_right
                    || !maybe_left
                        .borrow()
                        .is_bisect_left_of(&site)
                        .unwrap_or(false)
                {
                    break;
                }
            }
            let left = maybe_left.borrow().left_halfedge.as_ref().unwrap().clone();
            maybe_left = left;
        } else {
            // もしrightの可能性があるなら、leftに進む。
            loop {
                let left = maybe_left.borrow().left_halfedge.as_ref().unwrap().clone();
                maybe_left = left;

                // leftが本当にsiteの左なら、ループを止めてleftのrightを一番近いleftとしてみなす。
                // （もうやっているため、上のように別途指定しなくてもいい）
                if maybe_left.as_ptr() == p_most_left
                    || maybe_left
                        .borrow()
                        .is_bisect_left_of(&site)
                        .unwrap_or(false)
                {
                    break;
                }
            }
        }

        Some(maybe_left)
    }

    pub fn remove(&mut self, he: HalfEdgeRcCell) {
        // If left_halfedge is exist, 
        // we need to re-chain lhe->rhe = he.rhe.
        let mut he = he.borrow_mut();

        let rhe_cloned = he.right_halfedge.clone();
        let lhe_cloned = he.left_halfedge.clone();
        if let Some(lhe) = he.left_halfedge.as_mut() {
            lhe.borrow_mut().right_halfedge = rhe_cloned;
        }
        if let Some(rhe) = he.right_halfedge.as_mut() {
            rhe.borrow_mut().left_halfedge = lhe_cloned;
        }

        he.left_halfedge.take();
        he.right_halfedge.take();
    }

    pub fn visit_all<'a>(&'a self, func: &dyn Fn(std::cell::Ref<'_, HalfEdge>)) {
        let mut visited_id_set: HashSet<usize> = std::collections::HashSet::new();

        self.map.iter().for_each(|(_, ref_he)| {
            // Check this he is visited using id.
            let mut he = ref_he.clone();

            // If not visited yet, visit it.
            {
                let borrow_he = he.borrow();
                if !visited_id_set.contains(&borrow_he.id) {
                    visited_id_set.insert(borrow_he.id);
                    func(borrow_he);
                }
            }

            loop {
                // iterate chained right half-edges to end and visit if not visited.
                let next = match he.borrow().right_halfedge.as_ref() {
                    Some(next) => next.clone(),
                    None => break,
                };
                he = next;

                let borrow_he = he.borrow();
                if !visited_id_set.contains(&borrow_he.id) {
                    visited_id_set.insert(borrow_he.id);
                    func(borrow_he);
                }
            }
        })
    }
}

struct HalfEdgeVertexEvent {
    id: usize,
    vertex: FPoint2,
    y_star: f32,
    halfedge: HalfEdgeRcCell,
    /// `y_star`または`vertex[0] == x`が上がっていく。
    up_next: Option<HalfEdgeVertexEventRcCell>,
}

impl HalfEdgeVertexEvent {
    pub fn new(he: HalfEdgeRcCell, site: FPoint2, offset: f32) -> Self {
        Self {
            id: unsafe {
                let id = VE_ID_COUNTER;
                VE_ID_COUNTER += 1;
                id
            },
            vertex: site,
            y_star: site[1] + offset,
            halfedge: he,
            up_next: None,
        }
    }

    pub fn into_rccell(self) -> HalfEdgeVertexEventRcCell {
        Rc::new(RefCell::new(self))
    }
}

impl std::fmt::Debug for HalfEdgeVertexEvent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Just show left or right is exist.
        let halfedge_id = self.halfedge.borrow().id;
        let next_id = match self.up_next.as_ref() {
            Some(v) => Some(v.borrow().id),
            None => None,
        };

        f.debug_struct("VE")
            .field("ve_id", &self.id)
            .field("halfedge_id", &halfedge_id)
            .field("y_star", &self.y_star)
            .field("vertex", &self.vertex)
            .field("up_next", &next_id)
            .finish()
    }
}

struct HalfEdgeVertexEventMap {
    map: HalfEdgeVertexEventBTreeMap,
    bottom: HalfEdgeVertexEventRcCell,
}

impl HalfEdgeVertexEventMap {
    pub fn new() -> Self {
        Self {
            map: HalfEdgeVertexEventBTreeMap::new(),
            bottom: HalfEdgeVertexEvent::new(
                HalfEdge::default_rccell(),
                FPoint2::new(0f32, -5000f32),
                0f32,
            )
            .into_rccell(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    pub fn can_process_vertex_event_with(&self, site: FPoint2) -> bool {
        if self.is_empty() {
            return false;
        }

        match self.bottom.borrow().up_next.as_ref() {
            Some(min) => {
                let borrowed_min = min.borrow();
                let event_site = FPoint2::new(borrowed_min.vertex[0], borrowed_min.y_star);

                if site[1] < event_site[1] {
                    return false;
                } else if site[1] == event_site[1] && site[0] < event_site[0] {
                    return false;
                } else {
                    return true;
                }
            }
            None => false,
        }
    }

    pub fn try_extract_min(&mut self) -> Option<HalfEdgeVertexEventRcCell> {
        if self.is_empty() {
            None
        } else {
            // bottomとminのup_nextをつなげる。
            // min->up_nextは親を知らないので、bottomからつなげればいいだけ。
            let min = self.bottom.borrow_mut().up_next.as_ref().unwrap().clone();
            //bottom.up_next = min.borrow_mut().up_next.take();

            // minはまだ削除されていないので、btree-mapから削除する
            self.delete_halfedge(min.clone());

            Some(min)
        }
    }

    /// すでにmapの中に入ってる場合は考慮しない。
    pub fn insert_halfedge(&mut self, he: HalfEdgeRcCell, site: FPoint2, offset: f32) {
        let ve = HalfEdgeVertexEvent::new(he.clone(), site, offset);

        // `y`のkeyを計算する。
        let key = ve.y_star;

        // veよりyまたはxが大きいEventを探す。
        // もしマップに何もなければ、何もしない。(`None`のまま)
        let vecell = ve.into_rccell();
        if self.map.is_empty() {
            self.bottom.borrow_mut().up_next = Some(vecell.clone());
            he.borrow_mut().ve_ref = Some(vecell.clone());
            self.map.insert(OrderedFloat(key), vecell);

            return;
        }

        // 条件に合うlast(maybe_down)を探す。(veより手前下)
        let mut maybe_down = match self.map.range(..OrderedFloat(key)).next() {
            Some((_, v)) => v.clone(),
            None => self.bottom.clone(),
        };
        if maybe_down.borrow().up_next.is_some() {
            loop {
                let next = maybe_down.borrow().up_next.as_ref().unwrap().clone();
                let (next_y_star, vertex_x) = {
                    let bor_next = next.borrow();
                    (bor_next.y_star, bor_next.vertex[0])
                };
                if key > next_y_star || (key == next_y_star && site[0] > vertex_x) {
                    maybe_down = next;
                } else {
                    break;
                }
            }
        }

        // veのnextとmaybe_downのnextを再構築する。
        if let Some(next) = maybe_down.borrow().up_next.as_ref() {
            vecell.borrow_mut().up_next = Some(next.clone());
        }

        maybe_down.borrow_mut().up_next = Some(vecell.clone());
        he.borrow_mut().ve_ref = Some(vecell.clone());
        self.map.insert(OrderedFloat(key), vecell);
    }

    pub fn delete_halfedge(&mut self, he: HalfEdgeVertexEventRcCell) {
        let key = he.borrow().y_star;

        // 条件に合うlast(maybe_down)を探す。(veより手前下)
        let mut maybe_down = match self.map.range(..OrderedFloat(key)).next() {
            Some((_, v)) => v.clone(),
            None => self.bottom.clone(),
        };
        assert!(maybe_down.borrow().up_next.is_some());

        let p_he = he.as_ptr();
        loop {
            let p_down_next = maybe_down.borrow().up_next.as_ref().unwrap().as_ptr();
            if p_down_next != p_he {
                let next = maybe_down.borrow().up_next.as_ref().unwrap().clone();
                maybe_down = next;
            } else {
                break;
            }
        }

        self.map.remove(&OrderedFloat(key));
        match &he.borrow().up_next {
            Some(next) => maybe_down.borrow_mut().up_next = Some(next.clone()),
            None => maybe_down.borrow_mut().up_next = None,
        }
    }

    pub fn print_all(&self) {
        println!("{:?}", (-1, &self.bottom));
        self.map.iter().for_each(|(k, v)| {
            println!("{:?}", (k, v));
        });
    }
}

fn create_sorted_sites(delaunarys: &[FPoint2]) -> Option<Vec<FPoint2>> {
    // Make meta point list which contains triangle & edge meta information.
    let mut sites = delaunarys.iter().map(|p| *p).collect_vec();
    if sites.is_empty() {
        return None;
    }

    // 最初一番大きいyを持つポイントを最初に、そして大きいxを持つポイントを優先するようにソートする。
    sites.sort_by(|a, b| {
        use std::cmp::Ordering;
        match a[1].partial_cmp(&b[1]).unwrap() {
            Ordering::Equal => a[0].partial_cmp(&b[1]).unwrap(),
            y => y,
        }
    });
    Some(sites)
}

fn convert_to_voronoi(delaunarys: &[FPoint2]) -> Option<Vec<Edge>> {
    // Make meta point list which contains triangle & edge meta information.
    let mut site_queue: VecDeque<_> = {
        let sites = create_sorted_sites(delaunarys).unwrap();
        sites.iter().for_each(|t| println!("Sites : {:?}", t));
        sites.into_iter().collect()
    };

    // Make graph traversing meta_points.
    let bottom_site = site_queue.pop_front().unwrap();
    println!("bottom_site: {:?}", bottom_site);

    let mut halfedges = HalfEdgeMap::new();
    let mut vertex_events = HalfEdgeVertexEventMap::new();
    let voronoi_edges = Rc::new(RefCell::new(Vec::<Edge>::new()));

    // Check when inside of voronoi points.
    // sitesを全部通った後にも新しいbpが出来ることがある。
    let mut new_site: Option<FPoint2> = None;
    while !site_queue.is_empty() || !vertex_events.is_empty() {
        // Check we should process site event or voronoi vertex event.
        let is_site_event = if new_site.is_none() {
            true
        } else if site_queue.is_empty() {
            false
        } else {
            // Get new vertex event reference if available.
            let next_site = *site_queue.front().unwrap();
            !vertex_events.can_process_vertex_event_with(next_site)
        };

        //
        if is_site_event {
            let site = site_queue.pop_front().unwrap();
            new_site = Some(site.clone());

            println!("");
            println!("");
            println!("new_site : {:?}", new_site);

            // Get left half-edge and right half-edge boundary.
            // either given left or right may be end boundary half-edge.
            let mut l_boundary = halfedges.get_nearest_left_of(&site).unwrap().clone();
            let r_boundary = l_boundary.borrow().right_halfedge.as_ref().unwrap().clone();
            println!("l_bd : {:?}", l_boundary);
            println!("r_bd : {:?}", r_boundary);
            println!("");

            // Get left end from right boundary if exist, otherwise, just return bottom_site.
            let bot = l_boundary
                .borrow()
                .try_get_edge_end()
                .unwrap_or(bottom_site);
            let site_edge = Edge::new(bot, site);
            let mut l_halfedge = HalfEdge::from_edge(&site_edge, false).into_rccell();

            HalfEdge::chain_as_right(&mut l_boundary, l_halfedge.clone());
            println!("u l_bd : {:?}", l_boundary);
            println!("n l_he : {:?}", l_halfedge);
            println!("");

            // Check.. to create PQ for bisect.
            let is_leftbis_intersected = l_boundary
                .borrow()
                .try_get_bisect_intersected(&l_halfedge.borrow());
            if let Some(intersected) = is_leftbis_intersected {
                println!("lb-lh intersected : {:?}", intersected);
                {
                    let mut halfedge = l_boundary.borrow_mut();
                    if halfedge.ve_ref.is_some() {
                        let ve = halfedge.ve_ref.as_ref().unwrap().clone();
                        halfedge.ve_ref = None;
                        vertex_events.delete_halfedge(ve);
                    }
                }

                let dist = (intersected - site).magnitude();
                vertex_events.insert_halfedge(l_boundary.clone(), intersected, dist);
            };

            // Create reversed but dualed half edge.
            let ledge = l_halfedge.borrow().clone_edge().unwrap();
            let r_halfedge = HalfEdge::from_edge_container(ledge, true).into_rccell();
            HalfEdge::chain_as_right(&mut l_halfedge, r_halfedge.clone());
            //println!("u l_he : {:?}", l_halfedge);
            //println!("n r_he : {:?}", r_halfedge);
            //println!("");

            // Check intersect to create PQ for rev_bisect.
            let is_bisrev_intersected = r_halfedge
                .borrow()
                .try_get_bisect_intersected(&r_boundary.borrow());
            if let Some(intersected) = is_bisrev_intersected {
                println!("b-reb intersected : {:?}", intersected);

                let dist = (intersected - site).magnitude();
                vertex_events.insert_halfedge(r_halfedge.clone(), intersected, dist);
            };
        } else {
            // ここでvertex_eventが空っぽなのかを確認する必要はないかと。
            let (l_bnd, ve_point)  = {
                let lbnd_ve = vertex_events.try_extract_min().unwrap();
                let v = lbnd_ve.borrow().vertex;
                let l = lbnd_ve.borrow().halfedge.clone();
                (l, v)
            };
            //println!("");
            //println!("");
            //println!("new vertex event : {}, {:?}", ve_point, l_bnd);
            vertex_events.print_all();

            let (mut ll_bnd, mut r_bnd, bot) = {
                let borrowed_lbnd_he = l_bnd.borrow();
                (
                    borrowed_lbnd_he.left_halfedge.as_ref().unwrap().clone(),
                    borrowed_lbnd_he.right_halfedge.as_ref().unwrap().clone(),
                    borrowed_lbnd_he.try_get_edge_start().unwrap_or(bottom_site),
                )
            };
            let (rr_bnd, top) = {
                let borrowed_rbnd_he = r_bnd.borrow();
                (
                    borrowed_rbnd_he.right_halfedge.as_ref().unwrap().clone(),
                    borrowed_rbnd_he.try_get_edge_end().unwrap_or(bottom_site),
                )
            };
            //println!("lbnd : {:?}", l_bnd);
            //println!("llbnd : {:?}", ll_bnd);
            //println!("rbnd : {:?}", r_bnd);
            //println!("rrbnd : {:?}", rr_bnd);
            //println!("");

            // Update voronoi edge start - end points.
            // Check either l_bnd or r_bnd completes voronoi edge.
            {
                let ove = l_bnd.borrow_mut().update_voronoi_edge(ve_point, false);
                if let Some(ve) = ove {
                    println!("New Voronoi Edge (Closed) {:?}", ve);
                    voronoi_edges.borrow_mut().push(ve);
                }
            }
            {
                let ove = r_bnd.borrow_mut().update_voronoi_edge(ve_point, false);
                if let Some(ve) = ove {
                    println!("New Voronoi Edge (Closed) {:?}", ve);
                    voronoi_edges.borrow_mut().push(ve);
                }
            }

            // Remove half-edges and vertex-events.
            halfedges.remove(l_bnd);
            if let Some(ve) = r_bnd.borrow_mut().ve_ref.take() {
                vertex_events.delete_halfedge(ve);
            }
            halfedges.remove(r_bnd);

            // Create new half-edge, (bottom, top) (maybe inversely by the conditions)
            // Update voronoi edge into newly created he.
            let (bot, top, is_reversed) = if bot[1] > top[1] {
                (top, bot, true)
            } else {
                (bot, top, false)
            };
            let site_edge = Edge::new(bot, top);
            let mut new_he = HalfEdge::from_edge(&site_edge, is_reversed).into_rccell();

            HalfEdge::chain_as_right(&mut ll_bnd, new_he.clone());
            //println!("u ll_boundary : {:?}", ll_bnd);
            //println!("n new_halfedge : {:?}", new_he);
            //println!("");

            {
                let ove = new_he.borrow_mut().update_voronoi_edge(ve_point, true);
                if let Some(ve) = ove {
                    println!("New Voronoi Edge (Closed) {:?}", ve);
                    voronoi_edges.borrow_mut().push(ve);
                }
            }

            // Check and insert new vertex events if can.
            let is_llhe_nhe_intersected = ll_bnd
                .borrow()
                .try_get_bisect_intersected(&new_he.borrow());
            if let Some(intersected) = is_llhe_nhe_intersected {
                println!("ll-new intersected : {:?}", intersected);
                {
                    let mut halfedge = ll_bnd.borrow_mut();
                    if halfedge.ve_ref.is_some() {
                        let ve = halfedge.ve_ref.as_ref().unwrap().clone();
                        halfedge.ve_ref = None;
                        vertex_events.delete_halfedge(ve);
                    }
                }

                let dist = (intersected - bot).magnitude();
                vertex_events.insert_halfedge(ll_bnd.clone(), intersected, dist);
            };

            let is_nhe_rrhe_intersected = new_he
                .borrow()
                .try_get_bisect_intersected(&rr_bnd.borrow());
            if let Some(intersected) = is_nhe_rrhe_intersected {
                println!("new-rr intersected : {:?}", intersected);

                let dist = (intersected - bot).magnitude();
                vertex_events.insert_halfedge(new_he.clone(), intersected, dist);
            };
        }

        //halfedges.print_all();
        //vertex_events.print_all();

        // Check all half-edges are closed loop.
        // If not, assert because half-edge loop is invalid.
        halfedges.check_validation();
    }

    // 最後に完結していないVoronoi-edgeを吐き出す。
    // Half-edgeと１つのveの点を使って無限またはバウンダリーまで伸ばす。
    let opened_ves = voronoi_edges.clone();
    halfedges.visit_all(&|he| {
        let min_border = FPoint2::new(-100f32, -100f32);
        let max_border = FPoint2::new(100f32, 100f32);
        if he.edge.is_none() {
            return;
        }

        //println!("for he : {:?}", he);
        if let Some(ve) = he.try_get_voronoi_edge(min_border, max_border) {
            println!("New Voronoi Edge (Opened) {:?}", ve);
            opened_ves.borrow_mut().push(ve);
        }
    });

    let results = {
        let mut_ves = voronoi_edges.borrow_mut();
        mut_ves.iter().map(|&e| e).collect_vec()
    };
    Some(results)
}

fn main() {
    let points = [
        FPoint2::new(1f32, 0f32),
        FPoint2::new(-1f32, 0f32),
        FPoint2::new(1f32, 2f32),
        FPoint2::new(-1f32, 2f32),
        FPoint2::new(13.9f32, 6.76f32),
        FPoint2::new(12.7f32, 10.6f32),
        FPoint2::new(8.7f32, 7.7f32),
        FPoint2::new(7.1f32, 4.24f32),
        FPoint2::new(4.6f32, 11.44f32),
    ];

    // Voronoi edges using fortune's sweepline algorithm.
    let voronoi_edges = convert_to_voronoi(&points).unwrap();

    points.iter().for_each(|site| println!("Input site : {}", site));
    voronoi_edges.iter().enumerate().for_each(|(i, c)| println!("{:3}, Output edges : {:?}", i, c));
}
