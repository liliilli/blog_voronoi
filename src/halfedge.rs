use std::{
    cell::RefCell,
    collections::{BTreeMap, HashSet},
    rc::Rc,
};

use ordered_float::OrderedFloat;

use crate::{
    edge::{
        Direction, Edge, EdgeContainer, EdgeContainerRcCell, SiteEdge, SiteRcCell, VoronoiEdge,
        VoronoiEdgeType,
    },
    vertevent::HalfEdgeVertexEventRcCell,
    FPoint2,
};

static mut HE_ID_COUNTER: usize = 0;

#[derive(Default)]
pub struct HalfEdge {
    id: usize,
    left_halfedge: Option<HalfEdgeRcCell>,
    right_halfedge: Option<HalfEdgeRcCell>,
    pub edge: Option<EdgeContainerRcCell>,
    pub is_reversed_he: bool,
    pub ve_ref: Option<HalfEdgeVertexEventRcCell>,
}
///
pub type HalfEdgeRcCell = Rc<RefCell<HalfEdge>>;

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
            Some(ve) => Some(ve.borrow().id()),
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

    pub fn from_edge(edge: SiteEdge, is_reversed_he: bool) -> Self {
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

    pub fn id(&self) -> usize {
        self.id
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

    pub fn is_bisect_left_of(&self, point: FPoint2) -> Option<bool> {
        match &self.edge {
            Some(container) => {
                // 単純にbisectから左なのかを計算するのではなく、
                // reversedかではないかによって特殊に判定を行う必要がある。
                let c = container.borrow();
                let topsite = c.site_edge.end_point();

                let is_right_of_site = point[0] > topsite[0];
                if is_right_of_site && !self.is_reversed_he {
                    return Some(true);
                } else if !is_right_of_site && self.is_reversed_he {
                    return Some(false);
                }

                match c.is_bisect_left_of(&point) {
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

    pub fn try_get_edge_start(&self) -> Option<SiteRcCell> {
        match self.edge.as_ref() {
            Some(ec) => {
                let borrowed_ec = ec.borrow();
                Some(match self.is_reversed_he {
                    false => borrowed_ec.site_edge.clone_start(),
                    true => borrowed_ec.site_edge.clone_end(),
                })
            }
            None => None,
        }
    }

    pub fn try_get_edge_end(&self) -> Option<SiteRcCell> {
        match self.edge.as_ref() {
            Some(ec) => {
                let borrowed_ec = ec.borrow();
                Some(match self.is_reversed_he {
                    true => borrowed_ec.site_edge.clone_start(),
                    false => borrowed_ec.site_edge.clone_end(),
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
                        borrowed.site_edge.clone(),
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
                        borrowed.site_edge.clone(),
                        borrowed.bisector_pos,
                        borrowed.bisector_dir,
                    )
                }
                None => return None,
            }
        };

        // Check lhs_e and rhs_e has same parent.
        if lhs_e.end_point() == rhs_e.end_point() {
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
            let le = &lhs_e.end_point();
            let re = &rhs_e.end_point();

            if le[1] < re[1] || (le[1] == re[1] && le[0] < re[0]) {
                (self, &lhs_e)
            } else {
                (rhs, &rhs_e)
            }
        };

        // Right of site of e?
        if intersected[0] >= e.end_point()[0] {
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

    pub fn try_get_voronoi_edge(
        &self,
        min_boundary: FPoint2,
        max_boundary: FPoint2,
    ) -> Option<(Edge, VoronoiEdgeType)> {
        let (a, b, c) = {
            let ec = self.edge.as_ref().unwrap();
            let borrowed_ec = ec.borrow();
            //dbg!(&self);
            borrowed_ec.try_get_coefficients_of_bisect().unwrap()
        };
        if a.abs() < f32::EPSILON && b.abs() < f32::EPSILON {
            return None;
        }
        //dbg!((a, b, c));

        let (s1, s2) = {
            let ec = self.edge.as_ref().unwrap();
            let borrowed_ec = ec.borrow();
            match borrowed_ec.voronoi_edge {
                VoronoiEdge::InfFrom(s, e) => (s.clone(), e.clone()),
                _ => unreachable!(),
            }
        };
        let is_reversed = {
            match (s1, s2) {
                (None, None) => unreachable!(),
                (None, Some(_)) => Some(true),
                (Some(_), None) => Some(false),
                (Some(_), Some(_)) => None,
            }
        };

        // Get a line.
        if b.abs() > f32::EPSILON {
            let mut collision = None;
            let mut p1 = {
                let x = if let Some(p) = s1 {
                    if min_boundary[0] >= p[0] {
                        // left x collided
                        collision = Some(Direction::Left);
                        min_boundary[0]
                    } else {
                        p[0]
                    }
                } else {
                    // left x collided
                    collision = Some(Direction::Left);
                    min_boundary[0]
                };
                let y = ((a * x) + c) / b * -1f32;
                FPoint2::new(x, y)
            };
            let mut p2 = {
                let x = if let Some(p) = s2 {
                    if max_boundary[0] <= p[0] {
                        // Right x collided
                        collision = Some(Direction::Right);
                        max_boundary[0]
                    } else {
                        p[0]
                    }
                } else {
                    // Right x collided
                    collision = Some(Direction::Right);
                    max_boundary[0]
                };
                let y = ((a * x) + c) / b * -1f32;
                FPoint2::new(x, y)
            };
            //dbg!(p1, p2);

            // Check updated (p1, p2) is out of bound.
            let (miny, maxy) = (min_boundary[1], max_boundary[1]);
            if (p1[1] > maxy && p2[1] > maxy) || (p1[1] < miny && p2[1] < miny) {
                return None;
            }

            // Clip by y1->x1, and y2->x2.
            // a.is_subnormal()の場合、上の条件式で弾かれるのでaは有効な値があると仮定して勧めて良いかも。
            let p1y_clipped = p1[1].clamp(miny, maxy);
            if p1y_clipped != p1[1] {
                collision = match p1y_clipped == miny {
                    true => Some(Direction::Down),
                    false => Some(Direction::Up),
                };

                let x = ((b * p1y_clipped) + c) / a * -1f32;
                p1 = FPoint2::new(x, p1y_clipped);
            }

            let p2y_clipped = p2[1].clamp(miny, maxy);
            if p2y_clipped != p2[1] {
                collision = match p2y_clipped == miny {
                    true => Some(Direction::Down),
                    false => Some(Direction::Up),
                };

                let x = ((b * p2y_clipped) + c) / a * -1f32;
                p2 = FPoint2::new(x, p2y_clipped);
            }

            if let Some(col_dir) = collision {
                Some((
                    Edge::new(p1, p2),
                    VoronoiEdgeType::Opened((col_dir, is_reversed.unwrap_or(false))),
                ))
            } else {
                Some((Edge::new(p1, p2), VoronoiEdgeType::Closed))
            }
        } else {
            let mut collision = None;
            let mut p2 = {
                let y = if let Some(p) = s2 {
                    if min_boundary[1] >= p[1] {
                        // down y collided
                        collision = Some(Direction::Down);
                        min_boundary[1]
                    } else {
                        p[1]
                    }
                } else {
                    // down y collided
                    collision = Some(Direction::Down);
                    min_boundary[1]
                };
                let x = ((b * y) + c) / a * -1f32;
                FPoint2::new(x, y)
            };
            let mut p1 = {
                let y = if let Some(p) = s1 {
                    if max_boundary[1] <= p[1] {
                        // up y collided
                        collision = Some(Direction::Up);
                        max_boundary[1]
                    } else {
                        p[1]
                    }
                } else {
                    // up y collided
                    collision = Some(Direction::Up);
                    max_boundary[1]
                };
                let x = ((b * y) + c) / a * -1f32;
                FPoint2::new(x, y)
            };
            //dbg!(p1, p2);

            // Check updated (p1, p2) is out of bound.
            let (minx, maxx) = (min_boundary[0], max_boundary[0]);
            if (p1[0] > maxx && p2[0] > maxx) || (p1[0] < minx && p2[0] < minx) {
                return None;
            }

            // Clip by x1->y1, and x2->y2.
            // a.is_subnormal()の場合、上の条件式で弾かれるのでaは有効な値があると仮定して勧めて良いかも。
            let p1x_clipped = p1[0].clamp(minx, maxx);
            if p1x_clipped != p1[0] {
                collision = match p1x_clipped == minx {
                    true => Some(Direction::Left),
                    false => Some(Direction::Right),
                };

                let y = ((a * p1x_clipped) + c) / b * -1f32;
                p1 = FPoint2::new(p1x_clipped, y);
            }

            let p2x_clipped = p2[0].clamp(minx, maxx);
            if p2x_clipped != p2[0] {
                collision = match p2x_clipped == minx {
                    true => Some(Direction::Left),
                    false => Some(Direction::Right),
                };

                let y = ((a * p2x_clipped) + c) / b * -1f32;
                p2 = FPoint2::new(p2x_clipped, y);
            }

            if let Some(col_dir) = collision {
                Some((
                    Edge::new(p1, p2),
                    VoronoiEdgeType::Opened((col_dir, is_reversed.unwrap_or(false))),
                ))
            } else {
                Some((Edge::new(p1, p2), VoronoiEdgeType::Closed))
            }
        }
    }

    pub fn push_edge_to_sites(&mut self, ve: Edge, in_ve_type: VoronoiEdgeType) {
        match self.edge.as_mut() {
            Some(ec) => {
                let mut mut_ec = ec.borrow_mut();
                let ve_type = match in_ve_type {
                    // もしOpenedなら、HalfEdgeの方向によって逆にして入れる必要がある。
                    VoronoiEdgeType::Opened((dir, rev)) => match self.is_reversed_he {
                        true => VoronoiEdgeType::Opened((dir, !rev)),
                        false => VoronoiEdgeType::Opened((dir, !rev))
                    },
                    v => v,
                };

                mut_ec.site_edge.insert_edge(ve, ve_type);
            }
            _ => (),
        }
    }

    pub fn try_get_left_he(&self) -> Option<&HalfEdgeRcCell> {
        self.left_halfedge.as_ref()
    }

    pub fn try_get_right_he(&self) -> Option<&HalfEdgeRcCell> {
        self.right_halfedge.as_ref()
    }
}

#[derive(Debug)]
pub struct HalfEdgeMap {
    map: HalfEdgeBTreeMap,
    most_left: HalfEdgeRcCell,
    most_right: HalfEdgeRcCell,
}
/// KeyタイプはHalfEdgeの概略なx軸位置を示す。
/// 各Siteの一番近そうなHalfEdgeを取得するため。
pub type HalfEdgeBTreeMap = BTreeMap<OrderedFloat<f32>, HalfEdgeRcCell>;

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
            let site_edge = borrowed_cursor
                .edge
                .as_ref()
                .unwrap()
                .borrow()
                .site_edge
                .clone();

            match is_reversed_he {
                true => site_edge.end_point(),
                false => site_edge.start_point(),
            }
        };
        let goal_site = next_start_site;

        loop {
            if cursor.borrow().edge.is_none() {
                // Is most-right?
                break;
            }

            let (start, end) = {
                let borrowed_cursor = cursor.borrow();
                let is_reversed_he = borrowed_cursor.is_reversed_he;
                let site_edge = borrowed_cursor
                    .edge
                    .as_ref()
                    .unwrap()
                    .borrow()
                    .site_edge
                    .clone();

                match is_reversed_he {
                    false => (site_edge.start_point(), site_edge.end_point()),
                    true => (site_edge.end_point(), site_edge.start_point()),
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

    pub fn get_nearest_left_of<'a>(&'a self, site: FPoint2) -> Option<HalfEdgeRcCell> {
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
            && maybe_left.borrow().is_bisect_left_of(site).unwrap_or(true);

        if is_most_left || is_left_of_point {
            loop {
                let right = maybe_left.borrow().right_halfedge.as_ref().unwrap().clone();
                maybe_left = right;

                // rightが本当にsiteの右なら、ループを止めてrightのleftを一番近いleftとしてみなす。
                if maybe_left.as_ptr() == p_most_right
                    || !maybe_left.borrow().is_bisect_left_of(site).unwrap_or(false)
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
                    || maybe_left.borrow().is_bisect_left_of(site).unwrap_or(false)
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

    pub fn visit_all<'a>(&'a self, func: &dyn Fn(std::cell::RefMut<'_, HalfEdge>)) {
        let mut visited_id_set: HashSet<usize> = std::collections::HashSet::new();

        self.map.iter().for_each(|(_, ref_he)| {
            // Check this he is visited using id.
            let mut he = ref_he.clone();

            // If not visited yet, visit it.
            {
                let borrow_he = he.borrow_mut();
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

                let borrow_he = he.borrow_mut();
                if !visited_id_set.contains(&borrow_he.id) {
                    visited_id_set.insert(borrow_he.id);
                    func(borrow_he);
                }
            }
        })
    }
}
