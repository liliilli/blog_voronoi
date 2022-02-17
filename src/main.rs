#![feature(drain_filter)]

use std::{cell, collections::VecDeque, ops::Index};

use itertools::Itertools;
use nalgebra::{Matrix3, Point2, Point3, Vector2};
type FPoint2 = Point2<f32>;
type FPoint3 = Point3<f32>;
type FVector2 = Vector2<f32>;
type FMatrix3 = Matrix3<f32>;

#[derive(Clone, Debug)]
struct Triangle {
    pub points: [FPoint2; 3],
}

impl Triangle {
    pub fn new(first: FPoint2, second: FPoint2, third: FPoint2) -> Self {
        // まず３つの点が反時計歩行になるように整列する。
        // 奥の方に進むように…
        let first3 = FPoint3::new(first[0], 0f32, first[1]);
        let second3 = FPoint3::new(second[0], 0f32, second[1]);
        let third3 = FPoint3::new(third[0], 0f32, third[1]);

        let atob = second3 - first3;
        let btoc = third3 - second3;
        let normal = atob.cross(&btoc);

        let (a, b, c) = {
            if normal[1].is_sign_positive() {
                (second, first, third)
            } else {
                (first, second, third)
            }
        };

        Self { points: [a, b, c] }
    }

    pub fn get_circumcircle(&self) -> Circle {
        // Using cofactor expansion, get S_x, and S_y, and a.
        let a = self.a();
        let b = self.b();
        let c = self.c();
        let a_len_2 = (a - FPoint2::origin()).magnitude_squared();
        let b_len_2 = (b - FPoint2::origin()).magnitude_squared();
        let c_len_2 = (c - FPoint2::origin()).magnitude_squared();

        let s_x = FMatrix3::new(
            a_len_2, a[1], 1f32, b_len_2, b[1], 1f32, c_len_2, c[1], 1f32,
        )
        .determinant()
            * 0.5f32;

        let s_y = FMatrix3::new(
            a[0], a_len_2, 1f32, b[0], b_len_2, 1f32, c[0], c_len_2, 1f32,
        )
        .determinant()
            * 0.5f32;

        let adet =
            FMatrix3::new(a[0], a[1], 1f32, b[0], b[1], 1f32, c[0], c[1], 1f32).determinant();

        let bdet = FMatrix3::new(
            a[0], a[1], a_len_2, b[0], b[1], b_len_2, c[0], c[1], c_len_2,
        )
        .determinant();

        let s = FPoint2::new(s_x, s_y);
        let circumcenter = s / adet;
        let circumradius =
            ((bdet / adet) + ((s - FPoint2::origin()).magnitude_squared() / adet.powi(2))).sqrt();

        Circle {
            point: circumcenter,
            radius: circumradius,
        }
    }

    pub fn a(&self) -> FPoint2 {
        self.points[0]
    }
    pub fn b(&self) -> FPoint2 {
        self.points[1]
    }
    pub fn c(&self) -> FPoint2 {
        self.points[2]
    }

    pub fn get_edges(&self) -> [Edge; 3] {
        let a = self.a();
        let b = self.b();
        let c = self.c();

        [Edge::new(a, b), Edge::new(b, c), Edge::new(c, a)]
    }

    pub fn is_including_edge(&self, other: &Edge) -> bool {
        self.get_edges()
            .into_iter()
            .any(|target_edge| target_edge.is_nearly_same(other))
    }

    pub fn vertices(&self) -> &[FPoint2] {
        &self.points
    }

    fn sign(p1: &FPoint2, p2: &FPoint2, p3: &FPoint2) -> f32 {
        (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])
    }

    pub fn is_including_vertex(&self, other: &FPoint2) -> bool {
        let vertices = self.vertices();
        let (d1, d2, d3) = {
            let d1 = Self::sign(other, &vertices[1], &vertices[0]);
            let d2 = Self::sign(other, &vertices[2], &vertices[1]);
            let d3 = Self::sign(other, &vertices[0], &vertices[2]);
            (d1, d2, d3)
        };

        //println!("{:?} {d1}, {d2}, {d3}", &vertices);
        (d1 <= 0f32) && (d2 <= 0f32) && (d3 <= 0f32)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Circle {
    pub point: FPoint2,
    pub radius: f32,
}

impl Circle {
    pub fn is_including(&self, point: &FPoint2) -> bool {
        ((self.point - point).magnitude() - self.radius) <= 0f32
    }

    pub fn get_merged_circle(&self, other: &Self) -> Self {
        let to_self_offset = self.point - other.point;
        if !to_self_offset.magnitude_squared().is_normal() {
            // If point to point is subnormal, just update radius with maximum.
            Self {
                point: self.point,
                radius: self.radius.max(other.radius),
            }
        } else {
            // If direction can be calculated, get most far points between two circles.
            let to_self_dir = to_self_offset.normalize();
            let other_far_point = other.point - (to_self_dir * other.radius);
            let self_far_point = self.point + (to_self_dir * self.radius);

            let new_radius = (self_far_point - other_far_point).magnitude() * 0.5f32;
            let new_center = other_far_point + (to_self_dir * new_radius);
            Self {
                point: new_center,
                radius: new_radius,
            }
        }
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
}

pub fn is_nearly_same_fpoint2(lhs: FPoint2, rhs: FPoint2, epsilon: f32) -> bool {
    let x = lhs[0] - rhs[0];
    let y = lhs[1] - rhs[1];

    (x.abs() <= epsilon) && (y.abs() <= epsilon)
}

fn bowyer_watson(points: &[FPoint2]) -> Option<Vec<Triangle>> {
    if points.len() <= 2 {
        return None;
    }

    println!("Inputs : {:?}", points);

    // Add super-triangle to triangulation,
    // which is large enough to completely contain all the points in pointList.
    let super_triangle = {
        let mut min = points[0];
        let mut max = points[0];
        points.iter().for_each(|p| {
            min = FPoint2::new(min[0].min(p[0]), min[1].min(p[1]));
            max = FPoint2::new(max[0].max(p[0]), max[1].max(p[1]));
        });
        {
            let offset = FVector2::new(1e-1, 1e-1);
            min -= offset;
            max += offset;
        }

        let size = max - min;
        let top = {
            let offsetx = size[0] * 0.5f32;
            let offsety = offsetx * 60f32.to_radians().tan();
            max + FVector2::new(-offsetx, offsety)
        };
        let (left, right) = {
            let offsety = size[1];
            let offsetx = offsety / 60f32.to_radians().tan();
            let offset = FVector2::new(offsetx, 0f32);
            (min - offset, FPoint2::new(max[0], min[1]) + offset)
        };

        Triangle::new(top, left, right)
    };
    let mut triangulation: Vec<Triangle> = vec![super_triangle.clone()];

    for point in points.iter() {
        // Get bad triangles from triangulation.
        // Drained traingles from original list will be moved.
        let bad_triangles = triangulation
            .drain_filter(|triangle| triangle.get_circumcircle().is_including(point))
            .collect::<Vec<_>>();

        // Get polygons (valid edges from bad triangels) to check.
        let mut polygons = vec![];
        for i in 0..bad_triangles.len() {
            let source_bad = &bad_triangles[i];
            let source_edges = source_bad.get_edges();
            let mut addable_edges = source_edges
                .to_vec()
                .drain_filter(|edge| {
                    bad_triangles
                        .iter()
                        .enumerate()
                        .filter(|(ti, _)| *ti != i)
                        .any(|(_, target)| target.is_including_edge(edge))
                        == false
                })
                .collect::<Vec<_>>();
            polygons.append(&mut addable_edges);
        }

        // re-triangulate the polygonal hole.
        let mut new_triangles = polygons
            .into_iter()
            .filter(|edge| {
                let to_point = point - edge.start;
                let to_end = edge.end - edge.start;

                let a2 = to_end.magnitude_squared();
                let adotb = to_end.dot(&to_point);
                let reduct = to_point - (to_end * adotb * a2.recip());
                reduct.magnitude_squared().is_normal()
            })
            .map(|edge| Triangle::new(edge.start, edge.end, *point))
            .collect_vec();
        triangulation.append(&mut new_triangles);
    }

    let super_triangle_vertices = super_triangle.vertices();
    Some(
        triangulation
            .into_iter()
            .filter(|triangle| {
                super_triangle_vertices
                    .iter()
                    .any(|v| triangle.is_including_vertex(v))
                    == false
            })
            .collect_vec(),
    )
}

/// 方向性を持たないEdge。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct IndirectEdge {
    lhs: FPoint2,
    rhs: FPoint2,
}

impl IndirectEdge {
    pub fn from_edge(edge: &Edge) -> Self {
        Self::new(&edge.start, &edge.end)
    }

    pub fn new(one: &FPoint2, two: &FPoint2) -> Self {
        // one->twoが(-1, 0)方向へ向くようにして入れる。
        let dir = (two - one).normalize();
        let left_amount = FVector2::new(-1f32, 0f32).dot(&dir);
        if left_amount > 0f32 {
            return Self {
                lhs: *one,
                rhs: *two,
            };
        } else if left_amount < 0f32 {
            return Self {
                lhs: *two,
                rhs: *one,
            };
        }

        // もし90°であれば、下の方向を向くようにする。(0, -1)
        let down_amount = FVector2::new(0f32, -1f32).dot(&dir);
        if down_amount >= 0f32 {
            return Self {
                lhs: *one,
                rhs: *two,
            };
        } else {
            return Self {
                lhs: *two,
                rhs: *one,
            };
        }
    }

    pub fn get_center(&self) -> FPoint2 {
        self.lhs + ((self.rhs - self.lhs) * 0.5f32)
    }
}

impl std::hash::Hash for IndirectEdge {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let transmute_point = |p: &FPoint2| unsafe {
            let a = std::mem::transmute::<f32, u32>(p[0]);
            let b = std::mem::transmute::<f32, u32>(p[1]);
            (a, b)
        };
        let lhs = transmute_point(&self.lhs);
        let rhs = transmute_point(&self.rhs);
        lhs.hash(state);
        rhs.hash(state);
    }
}

impl std::cmp::Eq for IndirectEdge {}

#[derive(Clone, Debug, Eq, PartialEq, Default)]
struct EdgeRef {
    pub refs: Vec<(usize, usize)>,
}

impl EdgeRef {
    pub fn add(&mut self, triangle_i: usize, edge_i: usize) {
        self.refs.push((triangle_i, edge_i));
    }

    pub fn is_outside_of_diagram(&self) -> bool {
        self.refs.len() == 1
    }
}

#[derive(Debug, Clone)]
enum MedialAxis {
    InEdge(IndirectEdge),
    OutRay((FPoint2, FVector2)),
}

impl MedialAxis {
    pub fn is_inside(&self) -> bool {
        match self {
            MedialAxis::InEdge(_) => true,
            MedialAxis::OutRay(_) => false,
        }
    }

    pub fn as_inside(&self) -> Option<&IndirectEdge> {
        match self {
            MedialAxis::InEdge(ie) => Some(ie),
            MedialAxis::OutRay(_) => None,
        }
    }

    pub fn as_outside(&self) -> Option<(&FPoint2, &FVector2)> {
        match self {
            MedialAxis::InEdge(_) => None,
            MedialAxis::OutRay((p, d)) => Some((p, d)),
        }
    }
}

#[derive(Debug, Default)]
struct VoronoiCell {
    medial_axies: Vec<MedialAxis>,
}

impl VoronoiCell {
    pub fn is_medial_contained(&self, other: &MedialAxis) -> bool {
        if self.medial_axies.is_empty() {
            return false;
        }

        match other {
            MedialAxis::InEdge(ie) => self
                .medial_axies
                .iter()
                .filter(|ma| ma.is_inside())
                .find(|ma| ma.as_inside().unwrap() == ie)
                .is_some(),
            MedialAxis::OutRay((p, d)) => self
                .medial_axies
                .iter()
                .filter(|ma| !ma.is_inside())
                .find(|ma| {
                    let (tp, td) = ma.as_outside().unwrap();
                    tp == p && td == d
                })
                .is_some(),
        }
    }

    pub fn add_medial_axis(&mut self, medial_axis: MedialAxis) {
        self.medial_axies.push(medial_axis);
    }
}

#[derive(Debug)]
struct Temporary {
    pub ti: usize,
    pub indirect_edges: [IndirectEdge; 3],
    pub circumcircle_point: FPoint2,
}

fn convert_to_voronoi(delaunarys: &[Triangle]) -> Option<Vec<VoronoiCell>> {
    // まず三角形のそれぞれのEdgeがどの三角形につながっているかを情報として作る。
    // Edgeはdelaunarys[0],[1]を持っているなど…
    // そして三角形のそれぞれのEdgeは逆になっている可能性があるので、
    // それも踏まえて判別すべき。
    let mut edge_refs = std::collections::HashMap::new();
    delaunarys.iter().enumerate().for_each(|(ti, t)| {
        t.get_edges().into_iter().enumerate().for_each(|(ei, e)| {
            edge_refs
                .entry(IndirectEdge::from_edge(&e))
                .or_insert(EdgeRef::default())
                .add(ti, ei);
        })
    });

    // Make meta point list which contains triangle & edge meta information.
    let meta_points = delaunarys
        .iter()
        .enumerate()
        .map(|(ti, t)| Temporary {
            ti,
            indirect_edges: {
                let v = t
                    .get_edges()
                    .iter()
                    .map(|e| IndirectEdge::from_edge(&e))
                    .collect_vec();
                [v[0], v[1], v[2]]
            },
            circumcircle_point: t.get_circumcircle().point,
        })
        .collect_vec();
    if meta_points.is_empty() {
        return None;
    }
    meta_points.iter().for_each(|t| println!("Output temporary : {:?}", t));

    // Make graph traversing meta_points.
    let mut cells = Vec::<VoronoiCell>::new();
    for cursor in &meta_points {
        // 最初はどこが隣接しているか、どこが外枠なのかを調べる。
        let (outside, inside): (Vec<_>, Vec<_>) = cursor
            .indirect_edges
            .iter()
            .enumerate()
            .partition(|(_, ie)| edge_refs.get(ie).unwrap().is_outside_of_diagram());
        println!("Outside : {:?} || Inside : {:?}", outside, inside);

        let start = cursor.circumcircle_point;
        let is_start_outside = delaunarys[cursor.ti].is_including_vertex(&start);
        // voronoi_cellの位置をIndexとして保存する。
        let mut nearest_cells: [Option<usize>; 3] = [None; 3];

        // まずOutsideから。
        // Outsideからは自分のcircumcircle_pointとEdgeの中央を向くRayとして構成される。
        // circumcircle_pointから生成されるボロノイ図の部屋は2D上なら3つしかまで生成されない。
        outside.into_iter().for_each(|(_, ie)| {
            let medial_axis = {
                // ちなみにstartが三角形の中に入っていないこともあるので、
                // 確認して外側ならdirを逆にする。
                let to = ie.get_center();
                let ray_dir = match is_start_outside {
                    true => (to - start).normalize(),
                    false => (start - to).normalize(),
                };
                MedialAxis::OutRay((start, ray_dir))
            };

            // 同じmedial_axisを持つcellが存在しているかを見つける。（ここ最適化ポイント）
            // 外に進むmedial_axisは必ず2つ持つか全く無いか２択なので、判断しやすい。
            let exist_cells = cells
                .iter()
                .enumerate()
                .filter(|(_, c)| c.is_medial_contained(&medial_axis));
            let mut is_exist = false;
            for (ci, _) in exist_cells {
                is_exist = true;
                let slot = nearest_cells.iter_mut().find(|oci| oci.is_none()).unwrap();
                *slot = Some(ci);
            }
            if !is_exist 

                return;
            }

            // 見つからなかったら、2つ新しく作る、そしてmedial_axisを入れる。
            {
                let mut lhs = VoronoiCell::default();
                lhs.add_medial_axis(medial_axis.clone());
                cells.push(lhs);

                let slot = nearest_cells.iter_mut().find(|oci| oci.is_none()).unwrap();
                *slot = Some(cells.len());
            }
            {
                let mut rhs = VoronoiCell::default();
                rhs.add_medial_axis(medial_axis);
                cells.push(rhs);

                let slot = nearest_cells.iter_mut().find(|oci| oci.is_none()).unwrap();
                *slot = Some(cells.len());
            }
        });


        inside.into_iter().for_each(|(ei, ie)| {});
    }

    Some(cells)
}

fn main() {
    let triangles = bowyer_watson(&[
        FPoint2::new(1f32, 0f32),
        FPoint2::new(-1f32, 0f32),
        FPoint2::new(1f32, 2f32),
        FPoint2::new(-1f32, 2f32),
    ])
    .unwrap();
    triangles
        .iter()
        .for_each(|triangle| println!("Output triangle : {:?}", triangle));

    // Filter.
    triangles
        .iter()
        .map(|t| t.get_circumcircle())
        .unique_by(|c| unsafe {
            let p = &c.point;
            let ua = std::mem::transmute::<f32, u32>(p[0]);
            let ub = std::mem::transmute::<f32, u32>(p[1]);
            let ur = std::mem::transmute::<f32, u32>(c.radius);
            (ua, ub, ur)
        })
        .collect_vec()
        .iter()
        .for_each(|c| println!("Output circumcircle : {:?}", c));

    let cells = convert_to_voronoi(&triangles).unwrap();
    cells.iter().for_each(|c| println!("Output cell : {:?}", c));
}
