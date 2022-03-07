#![feature(drain_filter)]
mod edge;
mod halfedge;
mod vertevent;

use halfedge::{HalfEdge, HalfEdgeMap};
use vertevent::HalfEdgeVertexEventMap;

use std::{cell::RefCell, collections::VecDeque, rc::Rc};

use itertools::Itertools;
use nalgebra::Point2;
type FPoint2 = Point2<f32>;

use edge::{Edge, Site, SiteEdge, SiteRcCell, VoronoiEdgeType};

fn create_sorted_sites(delaunarys: &[FPoint2]) -> Option<Vec<SiteRcCell>> {
    // Make meta point list which contains triangle & edge meta information.
    let mut sites = delaunarys
        .iter()
        .map(|p| Site::from_point(*p).into_rccell())
        .collect_vec();
    if sites.is_empty() {
        return None;
    }

    // 最初一番大きいyを持つポイントを最初に、そして大きいxを持つポイントを優先するようにソートする。
    sites.sort_by(|a, b| {
        use std::cmp::Ordering;
        let ap = a.borrow().point;
        let bp = b.borrow().point;

        match ap[1].partial_cmp(&bp[1]).unwrap() {
            Ordering::Equal => ap[0].partial_cmp(&bp[1]).unwrap(),
            y => y,
        }
    });
    Some(sites)
}

fn convert_to_voronoi(delaunarys: &[FPoint2]) -> Option<(Vec<Edge>, Vec<SiteRcCell>)> {
    // Make meta point list which contains triangle & edge meta information.
    let mut sites = create_sorted_sites(delaunarys).unwrap();
    let mut site_queue: VecDeque<_> = {
        sites.iter().for_each(|t| println!("Sites : {:?}", t));
        sites.iter().map(|s| s.clone()).collect()
    };

    let bottom_site = site_queue.pop_front().unwrap();
    let mut halfedges = HalfEdgeMap::new();
    let mut vertex_events = HalfEdgeVertexEventMap::new();
    let voronoi_edges = Rc::new(RefCell::new(Vec::<Edge>::new()));

    // Check when inside of voronoi points.
    // sitesを全部通った後にも新しいbpが出来ることがある。
    let mut new_site: Option<SiteRcCell> = None;
    while !site_queue.is_empty() || !vertex_events.is_empty() {
        // Check we should process site event or voronoi vertex event.
        let is_site_event = if new_site.is_none() {
            true
        } else if site_queue.is_empty() {
            false
        } else {
            // Get new vertex event reference if available.
            let next_site_point = site_queue.front().unwrap().borrow().point;
            !vertex_events.can_process_vertex_event_with(next_site_point)
        };

        //
        if is_site_event {
            let site = site_queue.pop_front().unwrap();
            let site_point = site.borrow().point;
            new_site = Some(site.clone());
            //println!("");
            //println!("");
            //println!("new_site : {:?}", new_site);

            // Get left half-edge and right half-edge boundary.
            // either given left or right may be end boundary half-edge.
            let mut l_boundary = halfedges.get_nearest_left_of(site_point).unwrap().clone();
            let r_boundary = l_boundary.borrow().try_get_right_he().unwrap().clone();
            //dbg!(&new_site);

            // Get left end from right boundary if exist, otherwise, just return bottom_site.
            let bot = l_boundary
                .borrow()
                .try_get_edge_end()
                .unwrap_or(bottom_site.clone());
            let mut l_halfedge =
                HalfEdge::from_edge(SiteEdge::new(bot, site.clone()), false).into_rccell();

            HalfEdge::chain_as_right(&mut l_boundary, l_halfedge.clone());
            //dbg!(&l_boundary);
            //dbg!(&l_halfedge);
            //dbg!(&r_boundary);

            // Check.. to create PQ for bisect.
            let is_leftbis_intersected = l_boundary
                .borrow()
                .try_get_bisect_intersected(&l_halfedge.borrow());
            if let Some(intersected) = is_leftbis_intersected {
                //println!("lb-lh intersected : {:?}", intersected);
                {
                    let mut halfedge = l_boundary.borrow_mut();
                    if halfedge.ve_ref.is_some() {
                        let ve = halfedge.ve_ref.as_ref().unwrap().clone();
                        halfedge.ve_ref = None;
                        vertex_events.delete_halfedge(ve);
                    }
                }

                let dist = (intersected - site_point).magnitude();
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
                //println!("b-reb intersected : {:?}", intersected);

                let dist = (intersected - site_point).magnitude();
                vertex_events.insert_halfedge(r_halfedge.clone(), intersected, dist);
            };
        } else {
            // ここでvertex_eventが空っぽなのかを確認する必要はないかと。
            let (l_bnd, ve_point) = {
                let lbnd_ve = vertex_events.try_extract_min().unwrap();
                let v = lbnd_ve.borrow().vertex;
                let l = lbnd_ve.borrow().halfedge().clone();
                (l, v)
            };
            //println!("");
            //println!("");
            //println!("new vertex event : {}, {:?}", ve_point, l_bnd);
            //vertex_events.print_all();

            let (mut ll_bnd, mut r_bnd, bot) = {
                let borrowed_lbnd_he = l_bnd.borrow();
                (
                    borrowed_lbnd_he.try_get_left_he().unwrap().clone(),
                    borrowed_lbnd_he.try_get_right_he().unwrap().clone(),
                    borrowed_lbnd_he
                        .try_get_edge_start()
                        .unwrap_or(bottom_site.clone()),
                )
            };
            let (rr_bnd, top) = {
                let borrowed_rbnd_he = r_bnd.borrow();
                (
                    borrowed_rbnd_he.try_get_right_he().unwrap().clone(),
                    borrowed_rbnd_he
                        .try_get_edge_end()
                        .unwrap_or(bottom_site.clone()),
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
                    //println!("New Voronoi Edge (Closed) {:?} \n\tin HE {:?}", ve, l_bnd);
                    l_bnd
                        .borrow_mut()
                        .push_edge_to_sites(ve, VoronoiEdgeType::Closed);
                    voronoi_edges.borrow_mut().push(ve);
                }
            }
            {
                let ove = r_bnd.borrow_mut().update_voronoi_edge(ve_point, false);
                if let Some(ve) = ove {
                    //println!("New Voronoi Edge (Closed) {:?} \n\tin HE {:?}", ve, r_bnd);
                    r_bnd
                        .borrow_mut()
                        .push_edge_to_sites(ve, VoronoiEdgeType::Closed);
                    voronoi_edges.borrow_mut().push(ve);
                }
            }
            //dbg!(&l_bnd, &r_bnd);

            // Remove half-edges and vertex-events.
            halfedges.remove(l_bnd);
            if let Some(ve) = r_bnd.borrow_mut().ve_ref.take() {
                vertex_events.delete_halfedge(ve);
            }
            halfedges.remove(r_bnd);

            // Create new half-edge, (bottom, top) (maybe inversely by the conditions)
            // Update voronoi edge into newly created he.
            let (bot, top, is_reversed) = {
                let bot_point = bot.borrow().point;
                let top_point = top.borrow().point;
                if bot_point[1] > top_point[1] {
                    (top, bot, true)
                } else {
                    (bot, top, false)
                }
            };
            let bot_point = bot.borrow().point;
            let mut new_he =
                HalfEdge::from_edge(SiteEdge::new(bot, top), is_reversed).into_rccell();

            HalfEdge::chain_as_right(&mut ll_bnd, new_he.clone());
            //println!("u ll_boundary : {:?}", ll_bnd);
            //println!("n new_halfedge : {:?}", new_he);
            //println!("");

            {
                let ove = new_he.borrow_mut().update_voronoi_edge(ve_point, true);
                if let Some(ve) = ove {
                    //println!("New Voronoi Edge (Closed) {:?} \n\tin HE {:?}", ve, new_he);
                    new_he
                        .borrow_mut()
                        .push_edge_to_sites(ve, VoronoiEdgeType::Closed);
                    voronoi_edges.borrow_mut().push(ve);
                }
            }

            // Check and insert new vertex events if can.
            let is_llhe_nhe_intersected =
                ll_bnd.borrow().try_get_bisect_intersected(&new_he.borrow());
            if let Some(intersected) = is_llhe_nhe_intersected {
                //println!("ll-new intersected : {:?}", intersected);
                {
                    let mut halfedge = ll_bnd.borrow_mut();
                    if halfedge.ve_ref.is_some() {
                        let ve = halfedge.ve_ref.as_ref().unwrap().clone();
                        halfedge.ve_ref = None;
                        vertex_events.delete_halfedge(ve);
                    }
                }

                let dist = (intersected - bot_point).magnitude();
                vertex_events.insert_halfedge(ll_bnd.clone(), intersected, dist);
            };

            let is_nhe_rrhe_intersected =
                new_he.borrow().try_get_bisect_intersected(&rr_bnd.borrow());
            if let Some(intersected) = is_nhe_rrhe_intersected {
                //println!("new-rr intersected : {:?}", intersected);

                let dist = (intersected - bot_point).magnitude();
                vertex_events.insert_halfedge(new_he.clone(), intersected, dist);
            };
        }

        // Check all half-edges are closed loop.
        // If not, assert because half-edge loop is invalid.
        halfedges.check_validation();
    }

    // 最後に完結していないVoronoi-edgeを吐き出す。
    // Half-edgeと１つのveの点を使って無限またはバウンダリーまで伸ばす。
    //halfedges.print_all();
    let min_border = FPoint2::new(-100f32, -100f32);
    let max_border = FPoint2::new(100f32, 100f32);
    let opened_ves = voronoi_edges.clone();
    halfedges.visit_all(&|mut he| {
        if he.edge.is_none() {
            return;
        }

        if let Some((ve, ve_type)) = he.try_get_voronoi_edge(min_border, max_border) {
            he.push_edge_to_sites(ve, ve_type);
            opened_ves.borrow_mut().push(ve);
        }
    });

    // 各開いたボロノイ部屋にBoxバウンダリーからの線を入れる。
    let border_corners = {
        let border_left_corner = min_border;
        let border_up_corner = FPoint2::new(min_border[0], max_border[1]);
        let border_right_corner = max_border;
        let border_down_corner = FPoint2::new(max_border[0], min_border[1]);

        [
            (border_left_corner, border_up_corner),
            (border_up_corner, border_right_corner),
            (border_right_corner, border_down_corner),
            (border_down_corner, border_left_corner),
        ]
    };
    sites.iter_mut().for_each(|s| {
        let bmut_s = s.borrow_mut();


    }); 

    let ves_without_boundary= {
        let mut_ves = voronoi_edges.borrow_mut();
        mut_ves.iter().map(|&e| e).collect_vec()
    };
    Some((ves_without_boundary, sites))
}

fn main() {
    let points = [
        //FPoint2::new(-1f32, 1f32),
        //FPoint2::new(1f32, 1f32),
        //FPoint2::new(0f32, -1f32),

        //FPoint2::new(2f32, 0f32), FPoint2::new(-2f32, 0f32),
        //FPoint2::new(2f32, 2f32), FPoint2::new(-2f32, 2f32),

        FPoint2::new(1f32, 0f32), FPoint2::new(-1f32, 0f32),
        FPoint2::new(1f32, 2f32), FPoint2::new(-1f32, 2f32),
        FPoint2::new(13.9f32, 6.76f32), FPoint2::new(12.7f32, 10.6f32),
        FPoint2::new(8.7f32, 7.7f32), FPoint2::new(7.1f32, 4.24f32),
        FPoint2::new(4.6f32, 11.44f32),
    ];

    // Voronoi edges using fortune's sweepline algorithm.
    let (voronoi_edges, sites) = convert_to_voronoi(&points).unwrap();
    points
        .iter()
        .for_each(|site| println!("Input Site : {}", site));
    println!("");

    voronoi_edges
        .iter()
        .enumerate()
        .for_each(|(i, c)| println!("{:3}, Output : {:?}", i, c));
    println!("");

    sites.iter().enumerate().for_each(|(i, s)| {
        let bs = s.borrow();
        println!("{:3}, Site : {:?}, Closed: {:?}", i, bs.point, bs.is_closed());
        bs.voronoi_edges
            .iter()
            .enumerate()
            .for_each(|(i, &v)| println!("\t{:3}, {:?}", i, v));
    });
}

#[cfg(test)]
mod tests {
    use crate::{convert_to_voronoi, FPoint2};

    fn test_func(points: &[FPoint2]) {
        // Voronoi edges using fortune's sweepline algorithm.
        let (voronoi_edges, sites) = convert_to_voronoi(&points).unwrap();
        points
            .iter()
            .for_each(|site| println!("Input Site : {}", site));
        println!("");

        voronoi_edges
            .iter()
            .enumerate()
            .for_each(|(i, c)| println!("{:3}, Output : {:?}", i, c));
        println!("");

        sites.iter().enumerate().for_each(|(i, s)| {
            let bs = s.borrow();
            println!("{:3}, Site : {:?}", i, bs.point);
            bs.voronoi_edges
                .iter()
                .enumerate()
                .for_each(|(i, &v)| println!("\t{:3}, {:?}", i, v));
        });
    }

    #[test]
    fn case1() {
        let points = [
            FPoint2::new(2f32, 0f32), FPoint2::new(-2f32, 0f32),
            FPoint2::new(2f32, 2f32), FPoint2::new(-2f32, 2f32),
        ];
        test_func(&points);
    }

    #[test]
    fn case2() {
        let points = [
            FPoint2::new(1f32, 0f32), FPoint2::new(-1f32, 0f32),
            FPoint2::new(1f32, 2f32), FPoint2::new(-1f32, 2f32),
            FPoint2::new(13.9f32, 6.76f32), FPoint2::new(12.7f32, 10.6f32),
            FPoint2::new(8.7f32, 7.7f32), FPoint2::new(7.1f32, 4.24f32),
            FPoint2::new(4.6f32, 11.44f32),
        ];
        test_func(&points);
    }

    #[test]
    fn case3() {
        let points = [
            FPoint2::new(1f32, 0f32), FPoint2::new(-1f32, 0f32),
            FPoint2::new(1f32, 1f32), FPoint2::new(-1f32, 1f32),
            FPoint2::new(1f32, 2f32), FPoint2::new(-1f32, 2f32),
            FPoint2::new(1f32, 3f32), FPoint2::new(-1f32, 3f32),
            FPoint2::new(1f32, 4f32), FPoint2::new(-1f32, 4f32),
            FPoint2::new(1f32, 5f32), FPoint2::new(-1f32, 5f32),
            FPoint2::new(1f32, 6f32), FPoint2::new(-1f32, 6f32),
            FPoint2::new(1f32, 7f32), FPoint2::new(-1f32, 7f32),
            FPoint2::new(1f32, 8f32), FPoint2::new(-1f32, 8f32),
            FPoint2::new(1f32, 9f32), FPoint2::new(-1f32, 9f32),
            FPoint2::new(1f32, 10f32), FPoint2::new(-1f32, 10f32),
            FPoint2::new(1f32, 11f32), FPoint2::new(-1f32, 11f32),
            FPoint2::new(1f32, 12f32), FPoint2::new(-1f32, 12f32),
            FPoint2::new(1f32, 13f32), FPoint2::new(-1f32, 13f32),
            FPoint2::new(1f32, 14f32), FPoint2::new(-1f32, 14f32),
            FPoint2::new(1f32, 15f32), FPoint2::new(-1f32, 15f32),
            FPoint2::new(1f32, 16f32), FPoint2::new(-1f32, 16f32),
        ];
        test_func(&points);
    }

    #[test]
    fn case4() {
        let points = [
            FPoint2::new(1f32, 0f32), FPoint2::new(-1f32, 0f32),
            FPoint2::new(2f32, 1f32), FPoint2::new(-2f32, 1f32),
            FPoint2::new(3f32, 2f32), FPoint2::new(-3f32, 2f32),
            FPoint2::new(4f32, 3f32), FPoint2::new(-4f32, 3f32),
            FPoint2::new(3f32, 4f32), FPoint2::new(-3f32, 4f32),
            FPoint2::new(2f32, 5f32), FPoint2::new(-2f32, 5f32),
            FPoint2::new(1f32, 6f32), FPoint2::new(-1f32, 6f32),
            FPoint2::new(1f32, 7f32), FPoint2::new(-1f32, 7f32),
            FPoint2::new(1f32, 8f32), FPoint2::new(-1f32, 8f32),
            FPoint2::new(1f32, 9f32), FPoint2::new(-1f32, 9f32),
            FPoint2::new(2f32, 10f32), FPoint2::new(-2f32, 10f32),
            FPoint2::new(3f32, 11f32), FPoint2::new(-3f32, 11f32),
            FPoint2::new(4f32, 12f32), FPoint2::new(-4f32, 12f32),
            FPoint2::new(7f32, 13f32), FPoint2::new(-7f32, 13f32),
            FPoint2::new(6f32, 14f32), FPoint2::new(-6f32, 14f32),
            FPoint2::new(5f32, 15f32), FPoint2::new(-5f32, 15f32),
            FPoint2::new(4f32, 16f32), FPoint2::new(-4f32, 16f32),
        ];
        test_func(&points);
    }

    #[test]
    fn case5() {
        let points = [
            FPoint2::new(1f32, 0f32), FPoint2::new(-1f32, 0f32),
            FPoint2::new(2f32, 1f32), FPoint2::new(-2f32, 1f32),
            FPoint2::new(3f32, 2f32), FPoint2::new(-3f32, 2f32),
            FPoint2::new(4f32, 3f32), FPoint2::new(-4f32, 3f32),
            FPoint2::new(3f32, 4f32), FPoint2::new(-3f32, 4f32),
            FPoint2::new(2f32, 5f32), FPoint2::new(-2f32, 5f32),
            FPoint2::new(1f32, 6f32), FPoint2::new(-1f32, 6f32),
            FPoint2::new(1f32, 7f32), FPoint2::new(-1f32, 7f32),
            FPoint2::new(1f32, 8f32), FPoint2::new(-1f32, 8f32),
            FPoint2::new(1f32, 9f32), FPoint2::new(-1f32, 9f32),
            FPoint2::new(2f32, 10f32), FPoint2::new(-2f32, 10f32),
            FPoint2::new(3f32, 11f32), FPoint2::new(-3f32, 11f32),
            FPoint2::new(4f32, 12f32), FPoint2::new(-4f32, 12f32),
            FPoint2::new(7f32, 13f32), FPoint2::new(-7f32, 13f32),
            FPoint2::new(6f32, 14f32), FPoint2::new(-6f32, 14f32),
            FPoint2::new(5f32, 15f32), FPoint2::new(-5f32, 15f32),
            FPoint2::new(4f32, 16f32), FPoint2::new(-4f32, 16f32),
        ];
        test_func(&points);
    }

    #[test]
    fn case6() {
        let points = [
            FPoint2::new(1f32, 0f32), FPoint2::new(1f32, 2f32),
            FPoint2::new(2f32, 0f32), FPoint2::new(2f32, 2f32),
            FPoint2::new(3f32, 0f32), FPoint2::new(3f32, 2f32),
            FPoint2::new(4f32, 0f32), FPoint2::new(4f32, 2f32),
        ];
        test_func(&points);
    }
}
