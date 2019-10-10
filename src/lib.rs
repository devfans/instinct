use nalgebra::{self as na, *};
use std::cmp::{PartialOrd, Ordering};

pub trait InstinctUtils {
    // fn avg(&self, other: &Self) -> Self;
    fn instinct_eq(&self, other: &Self) -> bool;
    fn instinct_gt(&self, other: &Self) -> bool;
    fn instinct_ge(&self, other: &Self) -> bool;
    fn instinct_lt(&self, other: &Self) -> bool;
    fn instinct_le(&self, other: &Self) -> bool;

    fn instinct_zero(&self) -> bool;
    fn instinct_ord(&self, other: &Self) -> Ordering;
    fn instinct_delta_ord(&self) -> Ordering;
    fn instinct_ndelta_ord(&self) -> Ordering;
}

impl<N: RealField> InstinctUtils for N {
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
/// with edge cases:
/// two same points will give Equal not None
/// two same lines will give Equal not None
/// two same planes will give Equal not None
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
        {
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
            if $a.x != $b.x {
                list.sort_unstable_by(|a, b| a.x.partial_cmp(&b.x).unwrap());
            } else if $a.y != $b.y {
                list.sort_unstable_by(|a, b| a.y.partial_cmp(&b.y).unwrap());
            } else if $a.z != $b.z {
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
        {
            if same_point!($p, $a) {
                return false;
            }

            let line1 = ($b.coords - $a.coords).normalize();
            let line2 = ($c.coords - $a.coords).normalize();
            let target = ($p.coords - $a.coords).normalize();

            let norm = line1.cross(&line2);
            norm.dot(&target).instinct_zero()
        }
    };
    ($p: expr, $plane: expr) => {
        point_in_plane!($p, $plane.0, $plane.1, $plane.2)
    }
}

/// Sort point and plane represented by three points in ndc
macro_rules! sort_point_and_plane {
    ($p: expr, $a: expr, $b: expr, $c: expr) => {
        {
            if same_point!($p, $a) {
                return Ordering::Equal;
            }

            let line1 = ($b.coords - $a.coords).normalize();
            let line2 = ($c.coords - $a.coords).normalize();
            let target = ($p.coords - $a.coords).normalize();

            let norm = line1.cross(&line2);
            let delta = norm.dot(&target);
            
            if norm.z.instinct_zero() {
                Ordering::Equal
            } else if norm.z.is_positive() {
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

//--------------------------------
/// Compare two Point3
impl<N: RealField> InstinctOrd<Self> for Point3<N> {
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
impl<N: RealField> InstinctOrd<(&Point3<N>, &Point3<N>)> for Point3<N> {
    fn instinct_cmp(&self, other: &(&Point3<N>, &Point3<N>)) -> Ordering {
        if other.0.z.instinct_eq(&other.1.z) {
            // verticle with z - axis
            sort_points!(self, &other.0)
        } else if other.0.x.instinct_eq(&other.1.x) && other.0.y.instinct_eq(&other.1.y) {
            // parallel with z - axis
            Ordering::Equal
        } else {
            let v = self.z;
            let min = other.0.z.min(other.1.z);
            let max = other.0.z.max(other.1.z);
            if v.instinct_ge(&max) {
                return v.instinct_ord(&max);
            }
            if v.instinct_le(&min) {
                return v.instinct_ord(&min);
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
    fn instinct_cmp_ext(&self, other: &(&Point3<N>, &Point3<N>)) -> Option<Ordering> {
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

impl<N: RealField> InstinctOrd<(Point3<N>, Point3<N>)> for Point3<N> {
    fn instinct_cmp(&self, other: &(Point3<N>, Point3<N>)) -> Ordering {
        self.instinct_cmp(&(&other.0, &other.1))
    }
    fn instinct_cmp_ext(&self, other: &(Point3<N>, Point3<N>)) -> Option<Ordering> {
        self.instinct_cmp_ext(&(&other.0, &other.1))
    }
}


//--------------------------------
/// Compare point and a plane linked by three points
impl<N: RealField> InstinctOrd<(&Point3<N>, &Point3<N>, &Point3<N>)> for Point3<N> {
    fn instinct_cmp(&self, other: &(&Point3<N>, &Point3<N>, &Point3<N>)) -> Ordering {
        let v = self.z;
        let min = other.0.z.min(other.1.z).min(other.2.z);
        let max = other.0.z.max(other.1.z).max(other.2.z);
        if v.instinct_ge(&max) {
            return v.instinct_ord(&max);
        }
        if v.instinct_le(&min) {
            return v.instinct_ord(&min);
        }

        if !is_plane!(other) {
            let line = sort_line_points!(other);
            return self.instinct_cmp(&(line[0], line[2]));
        }

        let line1 = (other.1.coords - other.0.coords).normalize();
        let line2 = (other.2.coords - other.0.coords).normalize();
        let target = (self.coords - other.0.coords).normalize();

        let norm = line1.cross(&line2);
        let delta = norm.dot(&target);
        
        if norm.z.instinct_zero() {
            Ordering::Equal
        } else if norm.z.is_positive() {
            delta.instinct_delta_ord()
        } else {
            delta.instinct_delta_ord().reverse()
        }

    }
    fn instinct_cmp_ext(&self, other: &(&Point3<N>, &Point3<N>, &Point3<N>)) -> Option<Ordering> {
        if !is_plane!(other) {
            let line = sort_line_points!(other);
            return self.instinct_cmp_ext(&(line[0], line[2]));
        }

        let line1 = (other.1.coords - other.0.coords).normalize();
        let line2 = (other.2.coords - other.0.coords).normalize();
        let target = (self.coords - other.0.coords).normalize();

        let norm = line1.cross(&line2);
        let delta = norm.dot(&target);
        if delta.instinct_zero() {
            // Point in plane
            None
        } else if norm.z.instinct_zero() {
            // Plane parallel with z-axis
            Some(Ordering::Equal)
        } else if norm.z.is_positive() {
            Some(delta.instinct_delta_ord())
        } else {
            Some(delta.instinct_delta_ord().reverse())
        }
    }
}

impl<N: RealField> InstinctOrd<(Point3<N>, Point3<N>, Point3<N>)> for Point3<N> {
    fn instinct_cmp(&self, other: &(Point3<N>, Point3<N>, Point3<N>)) -> Ordering {
        self.instinct_cmp(&(&other.0, &other.1, &other.2))
    }
    fn instinct_cmp_ext(&self, other: &(Point3<N>, Point3<N>, Point3<N>)) -> Option<Ordering> {
        self.instinct_cmp_ext(&(&other.0, &other.1, &other.2))
    }
}


macro_rules! reverse_cmp_ext {
    ($self: expr, $other: expr) => {
        match $other.instinct_cmp_ext($self) {
            Some(res) => Some(res.reverse()),
            None => None
        }
    }
}

//--------------------------------
/// Compare line linked by two points and a point
impl<N: RealField> InstinctOrd<Point3<N>> for (Point3<N>, Point3<N>) {
    fn instinct_cmp(&self, other: &Point3<N>) -> Ordering {
        other.instinct_cmp(self).reverse()
    }
    fn instinct_cmp_ext(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp_ext!(self, other)
    }
}

impl<N: RealField> InstinctOrd<Point3<N>> for (&Point3<N>, &Point3<N>) {
    fn instinct_cmp(&self, other: &Point3<N>) -> Ordering {
        other.instinct_cmp(self).reverse()
    }
    fn instinct_cmp_ext(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp_ext!(self, other)
    }
}

//--------------------------------
/// Compare two lines linked by two points for each
impl<N: RealField> InstinctOrd<(&Point3<N>, &Point3<N>)> for (&Point3<N>, &Point3<N>) {
    fn instinct_cmp(&self, other: &(&Point3<N>, &Point3<N>)) -> Ordering {
        let v_min = self.0.z.min(self.1.z);
        let v_max = self.0.z.max(self.1.z);
        let min = other.0.z.min(other.1.z);
        let max = other.0.z.max(other.1.z);
        if v_min.instinct_ge(&max) {
            return v_min.instinct_ord(&max);
        }
        if v_max.instinct_le(&min) {
            return v_max.instinct_ord(&min);
        }
        
        if same_point!(self.0, self.1) {
            return self.0.instinct_cmp(other)
        } else if same_point!(other.0, other.1) {
            return other.0.instinct_cmp(self).reverse();
        }

        if is_plane!(self.0, other.0, other.1) {
        } else {
            self.1.instinct_cmp(other)
        }
    }
    fn instinct_cmp_ext(&self, other: &(&Point3<N>, &Point3<N>)) -> Option<Ordering> {
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
        } else if self.0.z.instinct_eq(&self.1.z) && other.0.z.instinct_eq(&other.1.z) {
            // Share same xy plane
            self.0.instinct_cmp_ext(other.0)
        } else {
            None
        }
    }
}

impl<N: RealField> InstinctOrd<(&Point3<N>, &Point3<N>)> for (Point3<N>, Point3<N>) {
    fn instinct_cmp(&self, other: &(&Point3<N>, &Point3<N>)) -> Ordering {
        (&self.0, &self.1).instinct_cmp(other)
    }
    fn instinct_cmp_ext(&self, other: &(&Point3<N>, &Point3<N>)) -> Option<Ordering> {
        (&self.0, &self.1).instinct_cmp_ext(other)
    }
}

impl<N: RealField> InstinctOrd<(Point3<N>, Point3<N>)> for (Point3<N>, Point3<N>) {
    fn instinct_cmp(&self, other: &(Point3<N>, Point3<N>)) -> Ordering {
        (&self.0, &self.1).instinct_cmp(&(&other.0, &other.1))
    }
    fn instinct_cmp_ext(&self, other: &(Point3<N>, Point3<N>)) -> Option<Ordering> {
        (&self.0, &self.1).instinct_cmp_ext(&(&other.0, &other.1))
    }
}

impl<N: RealField> InstinctOrd<(Point3<N>, Point3<N>)> for (&Point3<N>, &Point3<N>) {
    fn instinct_cmp(&self, other: &(Point3<N>, Point3<N>)) -> Ordering {
        self.instinct_cmp(&(&other.0, &other.1))
    }
    fn instinct_cmp_ext(&self, other: &(Point3<N>, Point3<N>)) -> Option<Ordering> {
        self.instinct_cmp_ext(&(&other.0, &other.1))
    }
}

/*
//--------------------------------
/// Compare line linked by two points and a point
impl<N: RealField> InstinctOrd<Point3<N>> for (Point3<N>, Point3<N>, Point3<N>) {
    fn instinct_cmp(&self, other: &Point3<N>) -> Option<Ordering> {
    }
}

impl<N: RealField> InstinctOrd<Point3<N>> for (&Point3<N>, &Point3<N>, &Point3<N>) {
    fn instinct_cmp(&self, other: &Point3<N>) -> Option<Ordering> {
    }
}
*/


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cmp_points() {
        assert_eq!(Point3::new(1.,2.,3.).instinct_cmp(&Point3::new(4.,5.,3.)), Some(Ordering::Equal));
        assert_eq!(Point3::new(1.,2.,3.00000001).instinct_cmp(&Point3::new(4.,5.,3.)), Some(Ordering::Greater));
        assert_eq!(Point3::new(1.,2.,2.999999).instinct_cmp(&Point3::new(4.,5.,3.)), Some(Ordering::Less));
    }

    #[test]
    fn test_cmp_point_line() {
        let line = (Point3::new(1., 2., 3.), Point3::new(5., 8., 9.));
        let mut point = Point3::new(1., 2., 3.);
        assert_eq!(line.instinct_cmp(&point), Some(Ordering::Equal));
        point = Point3::new(3., 5., 6.);
        assert_eq!(line.instinct_cmp(&point), Some(Ordering::Equal));
        point = Point3::new(3., 5., 6.0000001);
        assert_eq!(line.instinct_cmp(&point), Some(Ordering::Greater));
        point = Point3::new(3., 5., 5.9999999);
        assert_eq!(line.instinct_cmp(&point), Some(Ordering::Less));
    }

    #[test]
    fn test_cmp_point_plane() {
        let plane = (Point3::new(1., 2., 3.), Point3::new(1.4, 4., 6.), Point3::new(0., 0., 0.));
        let mut point = Point3::new(0., 0., 0.);
        assert_eq!(plane.instinct_cmp(&point), Some(Ordering::Equal));
        point = Point3::new(0.8, 2., 3.);
        assert_eq!(plane.instinct_cmp(&point), Some(Ordering::Equal));
        point = Point3::new(0.8, 2., 2.999999);
        assert_eq!(plane.instinct_cmp(&point), Some(Ordering::Greater));
        point = Point3::new(0.8, 2., 3.000001);
        assert_eq!(plane.instinct_cmp(&point), Some(Ordering::Less));
    }

    #[test]
    fn test_cmp_lines() {
        let a = (Point3::new(1., 2., 3.), Point3::new(5., 8., 9.));
        let b = (Point3::new(1., 2., 3.), Point3::new(5., 8., 9.));
        assert_eq!(a.instinct_cmp(&b), Some(Ordering::Equal));
    }
}



