#![feature(drain_filter)]
use ordered_float::OrderedFloat;

use std::{cell::RefCell, collections::BTreeMap, rc::Rc};

use itertools::Itertools;
use nalgebra::{Matrix3, Point2, Point3, Rotation2, Vector2};
type FPoint2 = Point2<f32>;
type FPoint3 = Point3<f32>;
type FVector2 = Vector2<f32>;
type FMatrix3 = Matrix3<f32>;

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
            None => None
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

    pub fn is_left_of_bisect(&self, point: &FPoint2) -> Option<bool> {
        match &self.edge {
            Some(container) => {
                // 単純にbisectから左なのかを計算するのではなく、
                // reversedかではないかによって特殊に判定を行う必要がある。
                let c = container.borrow();
                let is_right_of_site = point[0] > c.site_edge.end[0];
                if is_right_of_site && self.is_reversed_he == false {
                    return Some(true);
                }
                else if is_right_of_site == false && self.is_reversed_he == true {
                    return Some(false);
                }

                c.is_left_of_bisect(point)
            },
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
                    (borrowed.site_edge, borrowed.bisector_pos, borrowed.bisector_dir)
                },
                None => return None,
            }
        };
        let (rhs_e, rbv, rbd) = {
            match rhs.edge.as_ref() {
                Some(ec) => {
                    let borrowed = ec.borrow();
                    (borrowed.site_edge, borrowed.bisector_pos, borrowed.bisector_dir)
                },
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
        if s.is_finite() == false || t.is_finite() == false {
            return None;
        }
        let intersected = lbv + (lbd * s);

        // もしEdgeのend点が逆なら、ターゲットを逆にする。
        let (el, e) = {
            let le = &lhs_e.end;
            let re = &rhs_e.end;

            if le[1] < re[1] || (le[1] == re[1] && le[0] < re[0]) {
                (self, &lhs_e)
            }
            else {
                (rhs, &rhs_e)
            }
        };

        // Right of site of e?
        if intersected[0] >= e.end[0] {
            if el.is_reversed_he == false {
                return None;
            }
        }
        else {
            if el.is_reversed_he == true {
                return None;
            }
        }

        Some(intersected)
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
        let mut most_right = HalfEdge::default_rccell();
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
            && maybe_left
                .borrow()
                .is_left_of_bisect(&site)
                .unwrap_or(false);

        if is_most_left || is_left_of_point {
            loop {
                let right = maybe_left.borrow().right_halfedge.as_ref().unwrap().clone();
                maybe_left = right;

                // rightが本当にsiteの右なら、ループを止めてrightのleftを一番近いleftとしてみなす。
                if maybe_left.as_ptr() == p_most_right
                    || maybe_left
                        .borrow()
                        .is_left_of_bisect(&site)
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
                    || !maybe_left.borrow().is_left_of_bisect(&site).unwrap_or(true)
                {
                    break;
                }
            }
        }

        Some(maybe_left)
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

#[derive(Debug, Default)]
struct VoronoiCell {}

impl VoronoiCell {}

fn convert_to_voronoi(delaunarys: &[FPoint2]) -> Option<Vec<VoronoiCell>> {
    // Make meta point list which contains triangle & edge meta information.
    let sites = create_sorted_sites(delaunarys).unwrap();
    sites.iter().for_each(|t| println!("Sites : {:?}", t));

    // Make graph traversing meta_points.
    let bottom_site = *sites.first().unwrap();
    println!("bottom_site: {:?}", bottom_site);

    let mut halfedges = HalfEdgeMap::new();
    let mut vertex_events = HalfEdgeVertexEventMap::new();

    // Check when inside of voronoi points.
    // sitesを全部通った後にも新しいbpが出来ることがある。
    for new_site in sites.iter().skip(1) {
        println!("");
        println!("");
        println!("new_site : {:?}", new_site);

        // Get left half-edge and right half-edge boundary.
        // either given left or right may be end boundary half-edge.
        let mut l_boundary = halfedges.get_nearest_left_of(new_site).unwrap().clone();
        let r_boundary = l_boundary
            .borrow()
            .right_halfedge
            .as_ref()
            .unwrap()
            .clone();
        println!("l_boundary : {:?}", l_boundary);
        println!("r_boundary : {:?}", r_boundary);
        println!("");

        // Get left end from right boundary if exist, otherwise, just return bottom_site.
        let bot = match l_boundary.borrow().try_get_edge_end() {
            Some(bottom) => bottom,
            None => bottom_site,
        };
        let site_edge = Edge::new(bot, *new_site);
        let mut l_halfedge = HalfEdge::from_edge(&site_edge, false).into_rccell();

        HalfEdge::chain_as_right(&mut l_boundary, l_halfedge.clone());
        println!("u l_boundary : {:?}", l_boundary);
        println!("n l_halfedge : {:?}", l_halfedge);
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

            let dist = (intersected - new_site).magnitude();
            vertex_events.insert_halfedge(l_boundary.clone(), intersected, dist);
        };

        // Create reversed but dualed half edge.
        let ledge = l_halfedge.borrow().clone_edge().unwrap();
        let r_halfedge = HalfEdge::from_edge_container(ledge, true).into_rccell();
        HalfEdge::chain_as_right(&mut l_halfedge, r_halfedge.clone());
        println!("u l_he : {:?}", l_halfedge);
        println!("n r_he : {:?}", r_halfedge);
        println!("");

        // Check intersect to create PQ for rev_bisect.
        let is_bisrev_intersected = r_halfedge
            .borrow()
            .try_get_bisect_intersected(&r_boundary.borrow());
        if let Some(intersected) = is_bisrev_intersected {
            println!("b-reb intersected : {:?}", intersected);

            let dist = (intersected - new_site).magnitude();
            vertex_events.insert_halfedge(r_halfedge.clone(), intersected, dist);
        };

        halfedges.print_all();
    }

    let mut cells = Vec::<VoronoiCell>::new();
    Some(cells)
}

fn main() {
    let points = [
        FPoint2::new(1f32, 0f32),
        FPoint2::new(-1f32, 0f32),
        FPoint2::new(1f32, 2f32),
        FPoint2::new(-1f32, 2f32),
        //FPoint2::new(13.9f32, 6.76f32),
        //FPoint2::new(12.7f32, 10.6f32),
        //FPoint2::new(8.7f32, 7.7f32),
        //FPoint2::new(7.1f32, 4.24f32),
        //FPoint2::new(4.6f32, 11.44f32),
    ];

    // Dual
    let cells = convert_to_voronoi(&points).unwrap();
    cells.iter().for_each(|c| println!("Output cell : {:?}", c));

    //new_halfedge_map().into_iter().for_each(|v| println!("{:?}", v));
}
