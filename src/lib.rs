use nalgebra::{self as na, *};
use std::cmp::{PartialOrd, Ordering};

pub type Line<N> = (Point3<N>, Point3<N>);
pub type LineRef<'a, N> = (&'a Point3<N>, &'a Point3<N>);
pub type Plane<N> = (Point3<N>, Point3<N>, Point3<N>);
pub type PlaneRef<'a, N> = (&'a Point3<N>, &'a Point3<N>, &'a Point3<N>);

pub trait InstinctUtils {
    // fn avg(&self, other: &Self) -> Self;
    fn instinct_eq(&self, other: &Self) -> bool;
    fn instinct_gt(&self, other: &Self) -> bool;
    fn instinct_ge(&self, other: &Self) -> bool;
    fn instinct_lt(&self, other: &Self) -> bool;
    fn instinct_le(&self, other: &Self) -> bool;

    fn instinct_zero(&self) -> bool;
    fn instinct_not_negative(&self) -> bool;
    fn instinct_not_positive(&self) -> bool;
    fn instinct_positive(&self) -> bool;
    fn instinct_negative(&self) -> bool;
    fn instinct_ord(&self, other: &Self) -> Ordering;
    fn instinct_delta_ord(&self) -> Ordering;
    fn instinct_ndelta_ord(&self) -> Ordering;
}

impl<'a, N: RealField> InstinctUtils for N {
    fn instinct_eq(&self, other: &Self) -> bool {
        (*self - *other).abs() <= N::default_max_relative()
    }

    fn instinct_gt(&self, other: &Self) -> bool {
        (!self.instinct_eq(other)) && self > other
    }

    fn instinct_lt(&self, other: &Self) -> bool {
        (!self.instinct_eq(other)) && self < other
    }

    fn instinct_ge(&self, other: &Self) -> bool {
        self.instinct_eq(other) || self > other
    }

    fn instinct_le(&self, other: &Self) -> bool {
        self.instinct_eq(other) || self < other
    }

    fn instinct_zero(&self) -> bool {
        self.abs() <= N::default_max_relative()
    }

    fn instinct_not_negative(&self) -> bool {
        self.is_positive() || self.instinct_zero()
    }
    fn instinct_not_positive(&self) -> bool {
        self.is_negative() || self.instinct_zero()
    }

    fn instinct_positive(&self) -> bool {
        self.is_positive() && !self.instinct_zero()
    }
    fn instinct_negative(&self) -> bool {
        self.is_negative() && !self.instinct_zero()
    }
    fn instinct_ord(&self, other: &Self) -> Ordering {
        //TODO: NaN cmp can appear here? 
        if self.instinct_eq(other) {
            Ordering::Equal
        } else if self < other {
            Ordering::Greater
        } else {
            Ordering::Less
        }
    }

    fn instinct_delta_ord(&self) -> Ordering {
        if self.instinct_zero() {
            Ordering::Equal
        } else if self.is_positive() {
            Ordering::Greater
        } else {
            Ordering::Less
        }
    }

    fn instinct_ndelta_ord(&self) -> Ordering {
        self.instinct_delta_ord().reverse()
    }
}

/// instinct_cmp only compare points and finite(truncated) lines and polygons
/// Cross and Parallel cases would give a Equal
///
/// instinct_cmp_ext compare points and infinite lines and planes represented by polygons
/// returns None when cross happens and Equal in parallel situations
/// with edge cases below will gives Equal not None:
/// two same points
/// two same lines
/// two same planes
/// point in line
/// point in plane
/// line in plane
pub trait InstinctOrd<Rhs> {
    fn instinct_cmp(&self, other: &Rhs) -> Ordering;
    fn instinct_cmp_ext(&self, other: &Rhs) -> Option<Ordering>;
}

/// Check if same point
macro_rules! same_point {
    ($self: expr, $other: expr) => {
        $self.x.instinct_eq(&$other.x) &&
        $self.y.instinct_eq(&$other.y) &&
        $self.z.instinct_eq(&$other.z)
    };
    ($p0: expr, $p1: expr, $p2: expr) => {
        same_point!($p0, $p1) && same_point!($p1, $p2)
    }
}

/// Check if points share same  XY plane
macro_rules! same_z {
    ($a: expr, $b: expr, $c: expr) => {
        $a.z.instinct_eq(&$b.z) &&
        $a.z.instinct_eq(&$c.z) &&
        $b.z.instinct_eq(&$c.z)
    };
    ($p: expr) => {
        same_z!($p.0, $p.1, $p.2)
    }
}

/// Sort two points in ndc
macro_rules! sort_points {
    ($self: expr, $other: expr) => {
        if same_point!($self, $other) {
            Ordering::Equal
        } else {
            $self.z.instinct_ord(&$other.z)
        }
    }
}
/// Check if two line is parallel
macro_rules! is_lines_parallel {
    ($a: expr, $b: expr, $c: expr, $d: expr) => {
        {
            let v1 = ($b.coords - $a.coords).normalize();
            let v2 = ($d.coords - $c.coords).normalize();
            same_point!(&v1, &v2) || same_point!(&v1, -v2)
        }
    }; 
    ($a: expr, $b: expr, $c: expr) => {
        is_lines_parallel!($a.0, $a.1, $b, $c)
    };
    ($a: expr, $b: expr) => {
        is_lines_parallel!($a.0, $a.1, $b.0, $b.1)
    }
}

/// Check if point is in line
macro_rules! point_in_line {
    ($p: expr, $l0: expr, $l1: expr) => {
        if same_point!($p, $l0) || same_point!($p, $l1) {
            true
        } else {
            let v = $p.coords - $l0.coords;
            let v_line = $l1.coords - $l0.coords;
            let dot = v.dot(&v_line);
            (dot*dot).instinct_eq(&(v.norm_squared()*v_line.norm_squared()))
        }
    };
    ($p: expr, $l: expr) => {
        point_in_line!($p, $l.0, $l.1)
    }
}

/// Sort points in a line
macro_rules! sort_line_points {
    ($a: expr, $b: expr, $c: expr) => {
        {
            let mut list = [$a, $b, $c];
            if !$a.x.instinct_eq(&$b.x) {
                list.sort_unstable_by(|a, b| a.x.partial_cmp(&b.x).unwrap());
            } else if !$a.y.instinct_eq(&$b.y) {
                list.sort_unstable_by(|a, b| a.y.partial_cmp(&b.y).unwrap());
            } else if !$a.z.instinct_eq(&$b.z) {
                list.sort_unstable_by(|a, b| a.z.partial_cmp(&b.z).unwrap());
            }
            list
        }
    };
    ($p: expr) => {
        sort_line_points!($p.0, $p.1, $p.2)
    }
}

/// Check if a plane is represented by three points
macro_rules! is_plane {
    ($p0: expr, $p1: expr, $p2: expr) => {
        if same_point!($p0, $p1) {
            false
        } else {
            !point_in_line!($p2, $p0, $p1)
        }
    };
    ($p: expr) => {
        is_plane!($p.0, $p.1, $p.2)
    }
}

/// Check if a point is in the plane represented by three points
macro_rules! point_in_plane {
    ($p: expr, $a: expr, $b: expr, $c: expr) => {
        if same_point!($p, $a) {
            true
        } else {
            let line1 = ($b.coords - $a.coords).normalize();
            let line2 = ($c.coords - $a.coords).normalize();
            let target = ($p.coords - $a.coords).normalize();

            let normal = line1.cross(&line2);
            normal.dot(&target).instinct_zero()
        }
    };
    ($p: expr, $plane: expr) => {
        point_in_plane!($p, $plane.0, $plane.1, $plane.2)
    }
}

/// Sort point and plane represented by three points in ndc
macro_rules! sort_point_and_plane {
    ($p: expr, $a: expr, $b: expr, $c: expr) => {
        if same_point!($p, $a) {
            Ordering::Equal
        } else {
            let line1 = ($b.coords - $a.coords).normalize();
            let line2 = ($c.coords - $a.coords).normalize();
            let target = ($p.coords - $a.coords).normalize();

            let normal = line1.cross(&line2);
            let delta = normal.dot(&target);
            
            if normal.z.instinct_zero() {
                Ordering::Equal
            } else if normal.z.is_positive() {
                delta.instinct_delta_ord()
            } else {
                delta.instinct_delta_ord().reverse()
            }
        }
    };
    ($p: expr, $plane: expr) => {
        sort_point_and_plane!($p, $plane.0, $plane.1, $plane.2)
    }
}

/// Sort point and line in ndc
macro_rules! sort_point_and_line {
    ($p: expr, $l0: expr, $l1: expr) => {
        (($l1.z - $l0.z) * ($p.y - $l0.y) - ($l1.y - $l0.y) * ($p.z - $l0.z))
            .instinct_delta_ord()
    };
    ($p: expr, $l: expr) => {
        sort_point_and_line!($p, $l.0, $l.1)
    } 
}

/// Get normal vector of two lines
macro_rules! get_lines_normal {
    ($a: expr, $b: expr, $c: expr, $d: expr) => {
        {
            let d1 = $b.coords - $a.coords;
            let d2 = $d.coords - $c.coords;
            let (s, t) = closest_points_line_line_parameters($a, &d1, $c, &d2);
            let p1 = $a + d1*s;
            let p2 = $c + d2*t;
            p2.coords - p1.coords
        }
    };
    ($a: expr, $b: expr) => {
        get_lines_normal!($a.0, $a.1, $b.0, $b.1)
    }
}

/// Get normal vector of point and a plane
macro_rules! get_point_plane_normal {
    ($p: expr, $a: expr, $b: expr, $c: expr) => {
        {
            let d1 = $c.coords - $b.coords;
            let d2 = $b.coords - $a.coords;
            let (s, t) = closest_points_line_line_parameters($p, &d1, $a, &d2);
            let p1 = $p + d1*s;
            let p2 = $a + d2*t;
            p2.coords - p1.coords
        }
    };
    ($p: expr, $plane: expr) => {
        get_point_plane_normal!($p, $plane.0, $plane.1, $plane.2)
    }
}
 
/// Reverse Option<Ordering>
macro_rules! reverse_cmp_ext {
    ($self: expr, $other: expr) => {
        match $other.instinct_cmp_ext($self) {
            Some(res) => Some(res.reverse()),
            None => None
        }
    }
}

//--------------------------------
/// Compare two Point3
impl<'a, N: RealField> InstinctOrd<Self> for Point3<N> {
    fn instinct_cmp(&self, other: &Self) -> Ordering {
        // self.z.instinct_ord(&other.z)
        sort_points!(self, other)
    }
    fn instinct_cmp_ext(&self, other: &Self) -> Option<Ordering> {
        // Edge case: two same points will give Equal not None
        Some(sort_points!(self, other))
    }
}


//--------------------------------
/// Compare Point3 with a line linked by two points
///
/// It's comparing the point and the plane represented by the projection of the line on Z-Y plane
impl<'a, N: RealField> InstinctOrd<LineRef<'a, N>> for Point3<N> {
    fn instinct_cmp(&self, other: &LineRef<'a, N>) -> Ordering {
        if other.0.z.instinct_eq(&other.1.z) {
            sort_points!(self, &other.0)
        } else {
            let v = self.z;
            let min = other.0.z.min(other.1.z);
            let max = other.0.z.max(other.1.z);
            if v.instinct_ge(&max) {
                return Ordering::Greater;
            }
            if v.instinct_le(&min) {
                return Ordering::Less;
            }

            // Convert projection plane to YZ coords
            // let plane = Unit::new_unchecked(Vector2::new(other.1.z - other.0.z, other.1.y - other.0.y));
            // let target = Unit::new_unchecked(Vector2::new(self.z - other.0.z, self.y - other.0.y));
            // TODO: Do we really need to normalize the vector for more exact result?
            // let delta = plane.x * target.y - plane.y * target.x;
            // let delta = ($other.1.z - $other.0.z) * ($self.y - $other.0.y) - ($other.1.y - $other.0.y) * ($self.z - $other.0.z);
            // println!("line point {} {} {} {}", delta, delta.is_positive(), delta.is_negative(), delta.abs() <= N::zero());
            // delta.instinct_delta_ord()
            sort_point_and_line!(self, other)
        }
    }
    fn instinct_cmp_ext(&self, other: &LineRef<'a, N>) -> Option<Ordering> {
        // Same point
        if same_point!(&other.0, &other.1) {
           return Some(sort_points!(self, &other.0));
        }

        // Point in line
        if point_in_line!(self, other) {
            return None;
        }
        
        // Point not in line
        // let delta = ($other.1.z - $other.0.z) * ($self.y - $other.0.y) - ($other.1.y - $other.0.y) * ($self.z - $other.0.z);
        // println!("line point {} {} {} {}", delta, delta.is_positive(), delta.is_negative(), delta.abs() <= N::zero());
        // Some(delta.instinct_delta_ord())
        Some(sort_point_and_line!(self, other))

    }
}

impl<'a, N: RealField> InstinctOrd<Line<N>> for Point3<N> {
    fn instinct_cmp(&self, other: &Line<N>) -> Ordering {
        self.instinct_cmp(&(&other.0, &other.1))
    }
    fn instinct_cmp_ext(&self, other: &Line<N>) -> Option<Ordering> {
        self.instinct_cmp_ext(&(&other.0, &other.1))
    }
}


//--------------------------------
/// Compare point and a plane linked by three points
impl<'a, N: RealField> InstinctOrd<PlaneRef<'a, N>> for Point3<N> {
    fn instinct_cmp(&self, other: &PlaneRef<'a, N>) -> Ordering {
        if same_z!(other) {
            return self.z.instinct_ord(&other.0.z);
        }
        let v = self.z;
        let min = other.0.z.min(other.1.z).min(other.2.z);
        let max = other.0.z.max(other.1.z).max(other.2.z);
        if v.instinct_ge(&max) {
            return Ordering::Greater;
        }
        if v.instinct_le(&min) {
            return Ordering::Less;
        }

        if !is_plane!(other) {
            let line = sort_line_points!(other);
            return self.instinct_cmp(&(line[0], line[2]));
        }

        let line1 = (other.1.coords - other.0.coords).normalize();
        let line2 = (other.2.coords - other.0.coords).normalize();
        let target = (self.coords - other.0.coords).normalize();

        let normal = line1.cross(&line2);
        let delta = normal.dot(&target);
        
        if normal.z.instinct_zero() {
            Ordering::Equal
        } else if normal.z.is_positive() {
            delta.instinct_delta_ord()
        } else {
            delta.instinct_delta_ord().reverse()
        }

    }
    fn instinct_cmp_ext(&self, other: &PlaneRef<'a, N>) -> Option<Ordering> {
        if !is_plane!(other) {
            let line = sort_line_points!(other);
            return self.instinct_cmp_ext(&(line[0], line[2]));
        }

        let line1 = (other.1.coords - other.0.coords).normalize();
        let line2 = (other.2.coords - other.0.coords).normalize();
        let target = (self.coords - other.0.coords).normalize();

        let normal = line1.cross(&line2);
        let delta = normal.dot(&target);
        if delta.instinct_zero() {
            // Point in plane
            Some(Ordering::Equal)
        } else if normal.z.instinct_zero() {
            // Plane parallel with z-axis
            Some(Ordering::Equal)
        } else if normal.z.is_positive() {
            Some(delta.instinct_delta_ord())
        } else {
            Some(delta.instinct_delta_ord().reverse())
        }
    }
}

impl<N: RealField> InstinctOrd<Plane<N>> for Point3<N> {
    fn instinct_cmp(&self, other: &Plane<N>) -> Ordering {
        self.instinct_cmp(&(&other.0, &other.1, &other.2))
    }
    fn instinct_cmp_ext(&self, other: &Plane<N>) -> Option<Ordering> {
        self.instinct_cmp_ext(&(&other.0, &other.1, &other.2))
    }
}

impl<N: RealField> InstinctOrd<Point3<N>> for Plane<N> {
    fn instinct_cmp(&self, other: &Point3<N>) -> Ordering {
        other.instinct_cmp(&(&self.0, &self.1, &self.2)).reverse()
    }
    fn instinct_cmp_ext(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp_ext!(self, other)
    }
}

impl<'a, N: RealField> InstinctOrd<Point3<N>> for PlaneRef<'a, N> {
    fn instinct_cmp(&self, other: &Point3<N>) -> Ordering {
        other.instinct_cmp(self).reverse()
    }
    fn instinct_cmp_ext(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp_ext!(self, other)
    }
}


//--------------------------------
/// Compare line linked by two points and a point
impl<'a, N: RealField> InstinctOrd<Point3<N>> for Line<N> {
    fn instinct_cmp(&self, other: &Point3<N>) -> Ordering {
        other.instinct_cmp(self).reverse()
    }
    fn instinct_cmp_ext(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp_ext!(self, other)
    }
}

impl<'a, N: RealField> InstinctOrd<Point3<N>> for LineRef<'a, N> {
    fn instinct_cmp(&self, other: &Point3<N>) -> Ordering {
        other.instinct_cmp(self).reverse()
    }
    fn instinct_cmp_ext(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp_ext!(self, other)
    }
}

//--------------------------------
/// Compare two lines linked by two points for each
impl<'a, N: RealField> InstinctOrd<LineRef<'a, N>> for LineRef<'a, N> {
    fn instinct_cmp(&self, other: &LineRef<'a, N>) -> Ordering {
        if same_point!(self.0, self.1) {
            return self.0.instinct_cmp(other);
        } else if same_point!(other.0, other.1) {
            return other.0.instinct_cmp(self).reverse();
        }
        let v_min = self.0.z.min(self.1.z);
        let v_max = self.0.z.max(self.1.z);
        let min = other.0.z.min(other.1.z);
        let max = other.0.z.max(other.1.z);
        if v_min.instinct_eq(&max) {
            if v_max > max || min < v_min {
                return Ordering::Greater;
            }
        } else if v_min > max {
            return Ordering::Greater;
        }
        
        if v_max.instinct_eq(&min) {
            if v_min < max || max > v_max {
                return Ordering::Less;
            }
        } else if v_max < min {
            return Ordering::Less;
        }

        let o_0 = sort_point_and_line!(self.0, other);
        let o_1 = sort_point_and_line!(self.1, other);

        if point_in_line!(self.0, other) || o_0 == Ordering::Equal {
            self.1.instinct_cmp(other)
        } else if point_in_line!(self.1, other) || o_1 == Ordering::Equal {
            self.0.instinct_cmp(other)
        } else if o_0 == o_1 {
            o_0
        } else {
            get_lines_normal!(self, other).z.instinct_delta_ord().reverse()
        }
    }
    fn instinct_cmp_ext(&self, other: &LineRef<'a, N>) -> Option<Ordering> {
        if same_point!(self.0, self.1) {
            self.0.instinct_cmp_ext(other)
        } else if same_point!(other.0, other.1) {
            reverse_cmp_ext!(self, other.0)
        } else if is_lines_parallel!(self, other) {
            if point_in_line!(self.0, other) {
                Some(Ordering::Equal)
            } else {
                self.0.instinct_cmp_ext(other)
            }
        } else {
            let normal = get_lines_normal!(self, other);
            if normal.z.instinct_zero() {
                if normal.norm().instinct_zero() {
                    None
                } else {
                    Some(Ordering::Equal)
                }
            } else if normal.z.is_negative() {
                Some(Ordering::Greater)
            } else {
                Some(Ordering::Less)
            }
        }
    }
}

impl<'a, N: RealField> InstinctOrd<LineRef<'a, N>> for Line<N> {
    fn instinct_cmp(&self, other: &LineRef<'a, N>) -> Ordering {
        (&self.0, &self.1).instinct_cmp(other)
    }
    fn instinct_cmp_ext(&self, other: &LineRef<'a, N>) -> Option<Ordering> {
        (&self.0, &self.1).instinct_cmp_ext(other)
    }
}

impl<'a, N: RealField> InstinctOrd<Line<N>> for Line<N> {
    fn instinct_cmp(&self, other: &Line<N>) -> Ordering {
        (&self.0, &self.1).instinct_cmp(&(&other.0, &other.1))
    }
    fn instinct_cmp_ext(&self, other: &Line<N>) -> Option<Ordering> {
        (&self.0, &self.1).instinct_cmp_ext(&(&other.0, &other.1))
    }
}

impl<'a, N: RealField> InstinctOrd<Line<N>> for LineRef<'a, N> {
    fn instinct_cmp(&self, other: &Line<N>) -> Ordering {
        self.instinct_cmp(&(&other.0, &other.1))
    }
    fn instinct_cmp_ext(&self, other: &Line<N>) -> Option<Ordering> {
        self.instinct_cmp_ext(&(&other.0, &other.1))
    }
}

//--------------------------------
/// Compare a line and a plane
impl<'a, N: RealField> InstinctOrd<PlaneRef<'a, N>> for LineRef<'a, N> {
    fn instinct_cmp(&self, other: &PlaneRef<'a, N>) -> Ordering {
        if same_point!(self.0, self.1) {
            self.0.instinct_cmp(other)
        } else {
            let v_min = self.0.z.min(self.1.z);
            let v_max = self.0.z.max(self.1.z);
            let min = other.0.z.min(other.1.z).min(other.2.z);
            let max = other.0.z.max(other.1.z).max(other.2.z);
            if v_min.instinct_eq(&max) {
                if v_max > max || min < v_min {
                    return Ordering::Greater;
                }
            } else if v_min > max {
                return Ordering::Greater;
            }
            
            if v_max.instinct_eq(&min) {
                if v_min < max || max > v_max {
                    return Ordering::Less;
                }
            } else if v_max < min {
                return Ordering::Less;
            }

            if !is_plane!(other) {
                let line = sort_line_points!(other);
                self.instinct_cmp(&(line[0], line[2]))
            } else {
                let normal_0 = get_point_plane_normal!(self.0, other);
                let normal_1 = get_point_plane_normal!(self.1, other);
                if normal_0.dot(&normal_1).instinct_not_negative() {
                    if normal_0.z.abs() > normal_1.z.abs() {
                        normal_0.z.instinct_delta_ord().reverse()
                    } else {
                        normal_1.z.instinct_delta_ord().reverse()
                    }
                } else {
                    // TODO: Finish this
                    Ordering::Equal
                }
            }
        }
    }
    fn instinct_cmp_ext(&self, other: &PlaneRef<'a, N>) -> Option<Ordering> {
        if same_point!(self.0, self.1) {
            self.0.instinct_cmp_ext(other)
        } else if !is_plane!(other) {
            let line = sort_line_points!(other);
            self.instinct_cmp_ext(&(line[0], line[2]))
        } else {
            let c_0 = point_in_plane!(self.0, other);
            let c_1 = point_in_plane!(self.1, other);
            if c_0 {
                if c_1 {
                    Some(Ordering::Equal)
                } else {
                    None
                }
            } else if c_1 {
                None
            } else {
                let v_0 = other.0.coords - other.2.coords;
                let v_1 = other.1.coords - other.2.coords;
                let normal = v_0.cross(&v_1);
                let v = self.1.coords - self.0.coords;
                if normal.dot(&v).instinct_zero() {
                    Some(get_point_plane_normal!(self.0, other).z.instinct_delta_ord().reverse())
                } else {
                    None
                }
            }
        }
    }
}

impl<'a, N: RealField> InstinctOrd<PlaneRef<'a, N>> for Line<N> {
    fn instinct_cmp(&self, other: &PlaneRef<'a, N>) -> Ordering {
        (&self.0, &self.1).instinct_cmp(other)
    }
    fn instinct_cmp_ext(&self, other: &PlaneRef<'a, N>) -> Option<Ordering> {
        (&self.0, &self.1).instinct_cmp_ext(other)
    }
}
 
impl<'a, N: RealField> InstinctOrd<Plane<N>> for LineRef<'a, N> {
    fn instinct_cmp(&self, other: &Plane<N>) -> Ordering {
        self.instinct_cmp(&(&other.0, &other.1, &other.2))
    }
    fn instinct_cmp_ext(&self, other: &Plane<N>) -> Option<Ordering> {
        self.instinct_cmp_ext(&(&other.0, &other.1, &other.2))
    }
}

impl<N: RealField> InstinctOrd<Plane<N>> for Line<N> {
    fn instinct_cmp(&self, other: &Plane<N>) -> Ordering {
        (&self.0, &self.1).instinct_cmp(&(&other.0, &other.1))
    }
    fn instinct_cmp_ext(&self, other: &Plane<N>) -> Option<Ordering> {
        (&self.0, &self.1).instinct_cmp_ext(&(&other.0, &other.1))
    }
}

impl<'a, N: RealField> InstinctOrd<LineRef<'a, N>> for PlaneRef<'a, N> {
    fn instinct_cmp(&self, other: &LineRef<'a, N>) -> Ordering {
        other.instinct_cmp(self).reverse()
    }
    fn instinct_cmp_ext(&self, other: &LineRef<'a, N>) -> Option<Ordering> {
        reverse_cmp_ext!(self, other)
    }
}

impl<'a, N: RealField> InstinctOrd<Line<N>> for PlaneRef<'a, N> {
    fn instinct_cmp(&self, other: &Line<N>) -> Ordering {
        (&other.0, &other.1).instinct_cmp(self).reverse()
    }
    fn instinct_cmp_ext(&self, other: &Line<N>) -> Option<Ordering> {
        reverse_cmp_ext!(self, other)
    }
}

impl<'a, N: RealField> InstinctOrd<LineRef<'a, N>> for Plane<N> {
    fn instinct_cmp(&self, other: &LineRef<'a, N>) -> Ordering {
        other.instinct_cmp(&(&self.0, &self.1, &self.2)).reverse()
    }
    fn instinct_cmp_ext(&self, other: &LineRef<'a, N>) -> Option<Ordering> {
        reverse_cmp_ext!(self, other)
    }
}
impl<N: RealField> InstinctOrd<Line<N>> for Plane<N> {
    fn instinct_cmp(&self, other: &Line<N>) -> Ordering {
        (&other.0, &other.1).instinct_cmp(&(&self.0, &self.1, &self.2)).reverse()
    }
    fn instinct_cmp_ext(&self, other: &Line<N>) -> Option<Ordering> {
        reverse_cmp_ext!(self, other)
    }
}

//--------------------------------
/// Compare two planes
impl<'a, N: RealField> InstinctOrd<PlaneRef<'a, N>> for PlaneRef<'a, N> {
    fn instinct_cmp(&self, other: &PlaneRef<'a, N>) -> Ordering {
    }
    fn instinct_cmp_ext(&self, other: &PlaneRef<'a, N>) -> Option<Ordering> {
        if !is_plane!(self) {
            let line = sort_line_points!(self);
            (line[0], line[2]).instinct_cmp_ext(other)
        } else if !is_plane!(other) {
            let line = sort_line_points!(other);
            self.instinct_cmp_ext(&(line[0], line[2]))
        } else {
            let n_0 = get_point_plane_normal!(self.0, other);
            let n_1 = get_point_plane_normal!(self.1, other);
            let n_2 = get_point_plane_normal!(self.2, other);
            if same_point!(&n_0, &n_1, &n_2) {
                Some(n_0.z.instinct_delta_ord().reverse())
            } else {
                None
            }
        }
    }
}

impl<'a, N: RealField> InstinctOrd<Plane<N>> for PlaneRef<'a, N> {
    fn instinct_cmp(&self, other: &Plane<N>) -> Ordering {
        self.instinct_cmp(&(&other.0, &other.1, &other.2))
    }
    fn instinct_cmp_ext(&self, other: &Plane<N>) -> Option<Ordering> {
        self.instinct_cmp_ext(&(&other.0, &other.1, &other.2))
    }
}
 
impl<'a, N: RealField> InstinctOrd<PlaneRef<'a, N>> for Plane<N> {
    fn instinct_cmp(&self, other: &PlaneRef<'a, N>) -> Ordering {
        (&self.0, &self.1, &self.2).instinct_cmp(other)
    }
    fn instinct_cmp_ext(&self, other: &PlaneRef<'a, N>) -> Option<Ordering> {
        (&self.0, &self.1, &self.2).instinct_cmp_ext(other)
    }
}
 
impl<N: RealField> InstinctOrd<Plane<N>> for Plane<N> {
    fn instinct_cmp(&self, other: &Plane<N>) -> Ordering {
        (&self.0, &self.1, &self.2).instinct_cmp(&(&other.0, &other.1, &other.2))
    }
    fn instinct_cmp_ext(&self, other: &Plane<N>) -> Option<Ordering> {
        (&self.0, &self.1, &self.2).instinct_cmp_ext(&(&other.0, &other.1, &other.2))
    }
}


//--------------------Utils borrowed from ncollide

/// Closest points between two lines.
///
/// The result, say `res`, is such that the closest points between both lines are
/// `orig1 + dir1 * res.0` and `orig2 + dir2 * res.1`.
#[inline]
pub fn closest_points_line_line_parameters<N: RealField>(
    orig1: &Point3<N>,
    dir1: &Vector3<N>,
    orig2: &Point3<N>,
    dir2: &Vector3<N>,
) -> (N, N) {
    let res =
        closest_points_line_line_parameters_eps(orig1, dir1, orig2, dir2, N::default_epsilon());
    (res.0, res.1)
}

/// Closest points between two lines with a custom tolerance epsilon.
///
/// The result, say `res`, is such that the closest points between both lines are
/// `orig1 + dir1 * res.0` and `orig2 + dir2 * res.1`. If the lines are parallel
/// then `res.2` is set to `true` and the returned closest points are `orig1` and
/// its projection on the second line.
#[inline]
pub fn closest_points_line_line_parameters_eps<N: RealField>(
    orig1: &Point3<N>,
    dir1: &Vector3<N>,
    orig2: &Point3<N>,
    dir2: &Vector3<N>,
    eps: N,
) -> (N, N, bool) {
    // Inspired by RealField-time collision detection by Christer Ericson.
    let r = *orig1 - *orig2;

    let a = dir1.norm_squared();
    let e = dir2.norm_squared();
    let f = dir2.dot(&r);

    let _0: N = na::zero();
    let _1: N = na::one();

    if a <= eps && e <= eps {
        (_0, _0, false)
    } else if a <= eps {
        (_0, f / e, false)
    } else {
        let c = dir1.dot(&r);
        if e <= eps {
            (-c / a, _0, false)
        } else {
            let b = dir1.dot(dir2);
            let ae = a * e;
            let bb = b * b;
            let denom = ae - bb;

            // Use absolute and ulps error to test collinearity.
            // let parallel = denom <= eps || ulps_eq!(ae, bb);
            let parallel = denom <= eps;

            let s = if !parallel {
                (b * f - c * e) / denom
            } else {
                _0
            };

            (s, (b * s + f) / e, parallel)
        }
    }
}

// FIXME: can we re-used this for the segment/segment case?
/// Closest points between two segments.
#[inline]
pub fn closest_points_line_line<N: RealField>(
    orig1: &Point3<N>,
    dir1: &Vector3<N>,
    orig2: &Point3<N>,
    dir2: &Vector3<N>,
) -> Line<N> {
    let (s, t) = closest_points_line_line_parameters(orig1, dir1, orig2, dir2);
    (*orig1 + *dir1 * s, *orig2 + *dir2 * t)
}


//-----------------test cases-------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cmp_points() {
        assert_eq!(Point3::new(1.,2.,3.).instinct_cmp(&Point3::new(4.,5.,3.)), Ordering::Equal);
        assert_eq!(Point3::new(1.,2.,3.00000001).instinct_cmp(&Point3::new(4.,5.,3.)), Ordering::Greater);
        assert_eq!(Point3::new(1.,2.,2.999999).instinct_cmp(&Point3::new(4.,5.,3.)), Ordering::Less);
    }

    #[test]
    fn test_cmp_point_line() {
        let line = (Point3::new(1., 2., 3.), Point3::new(5., 8., 9.));
        let mut point = Point3::new(1., 2., 3.);
        assert_eq!(line.instinct_cmp(&point), Ordering::Equal);
        point = Point3::new(3., 5., 6.);
        assert_eq!(line.instinct_cmp(&point), Ordering::Equal);
        point = Point3::new(3., 5., 6.0000001);
        assert_eq!(line.instinct_cmp(&point), Ordering::Greater);
        point = Point3::new(3., 5., 5.9999999);
        assert_eq!(line.instinct_cmp(&point), Ordering::Less);
    }

    #[test]
    fn test_cmp_point_plane() {
        let plane = (Point3::new(1., 2., 3.), Point3::new(1.4, 4., 6.), Point3::new(0., 0., 0.));
        let mut point = Point3::new(0., 0., 0.);
        assert_eq!(plane.instinct_cmp(&point), Ordering::Equal);
        point = Point3::new(0.8, 2., 3.);
        assert_eq!(plane.instinct_cmp(&point), Ordering::Equal);
        point = Point3::new(0.8, 2., 2.999999);
        assert_eq!(plane.instinct_cmp(&point), Ordering::Greater);
        point = Point3::new(0.8, 2., 3.000001);
        assert_eq!(plane.instinct_cmp(&point), Ordering::Less);
    }

    #[test]
    fn test_cmp_lines() {
        let a = (Point3::new(1., 2., 3.), Point3::new(5., 8., 9.));
        let b = (Point3::new(1., 2., 3.), Point3::new(5., 8., 9.));
        assert_eq!(a.instinct_cmp(&b), Ordering::Equal);
    }

    #[test]
    fn test_lines_perpendicular_vector() {
        let p1 = Point3::new(-1., 0., -1.);
        let p2 = Point3::new(1., 32., 0.);
        let d1 = Vector3::new(-10., 0., -10.);
        let d2 = Vector3::new(3., 0., 0.);
        let (s, t) = closest_points_line_line_parameters(&p1, &d1, &p2, &d2);
        let t1 = p1 + d1*s;
        let t2 = p2 + d2*t;
        assert_eq!(t1, Point3::new(0., 0., 0.));
        assert_eq!(t2, Point3::new(0., 32., 0.));
    }
}



