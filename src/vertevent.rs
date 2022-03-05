use std::{cell::RefCell, collections::BTreeMap, rc::Rc};

use ordered_float::OrderedFloat;

use crate::{
    halfedge::{HalfEdge, HalfEdgeRcCell},
    FPoint2,
};

///
pub type HalfEdgeVertexEventRcCell = Rc<RefCell<HalfEdgeVertexEvent>>;
pub type HalfEdgeVertexEventBTreeMap = BTreeMap<OrderedFloat<f32>, HalfEdgeVertexEventRcCell>;

static mut VE_ID_COUNTER: usize = 0;

pub struct HalfEdgeVertexEvent {
    id: usize,
    pub vertex: FPoint2,
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

    pub fn id(&self) -> usize {
        self.id
    }

    pub fn halfedge(&self) -> &HalfEdgeRcCell {
        &self.halfedge
    }
}

impl std::fmt::Debug for HalfEdgeVertexEvent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Just show left or right is exist.
        let halfedge_id = self.halfedge.borrow().id();
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

pub struct HalfEdgeVertexEventMap {
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
