use nalgebra::{self as na, *};
use std::cmp::{PartialOrd, Ordering};

pub trait Utils {
    fn avg(&self, other: &Self) -> Self;
}

pub trait InstinctOrd<Rhs> {
    fn instinct_cmp(&self, other: &Rhs) -> Option<Ordering>;
}

//--------------------------------
/// Compare two Point3
impl<N: Scalar + PartialOrd> InstinctOrd<Self> for Point3<N> {
    fn instinct_cmp(&self, other: &Self) -> Option<Ordering> {
        self.z.partial_cmp(&other.z)
    }
}

//--------------------------------
/// Compare Point3 with a line linked by two points
///
/// It's comparing the point and the plane represented by the projection of the line on Z-Y plane
macro_rules! impl_point_cmp_line {
    ($self: ident, $other: ident) => {
        if $other.0.z == $other.1.z {
            $self.z.partial_cmp(&$other.0.z)
        } else if $other.0.x == $other.1.x && $other.0.y == $other.1.y {
            Some(Ordering::Equal)
        } else {
            let v = $self.z;
            let min = $other.0.z.min($other.1.z);
            let max = $other.0.z.max($other.1.z);
            if v >= max {
                return v.partial_cmp(&max);
            }
            if v <= min {
                return v.partial_cmp(&min);
            }

            // Convert projection plane to YZ coords
            // let plane = Unit::new_unchecked(Vector2::new(other.1.z - other.0.z, other.1.y - other.0.y));
            // let target = Unit::new_unchecked(Vector2::new(self.z - other.0.z, self.y - other.0.y));
            // TODO: Do we really need to convert to unit vector for more exact result?
            // let delta = plane.x * target.y - plane.y * target.x;
            let delta = ($other.1.z - $other.0.z) * ($self.y - $other.0.y) - ($other.1.y - $other.0.y) * ($self.z - $other.0.z);
            if delta.is_positive() {
                Some(Ordering::Greater)
            } else if delta.is_negative() {
                Some(Ordering::Less)
            } else {
                Some(Ordering::Equal)
            }
        }
    }
}
impl<N: RealField> InstinctOrd<(&Point3<N>, &Point3<N>)> for Point3<N> {
    fn instinct_cmp(&self, other: &(&Point3<N>, &Point3<N>)) -> Option<Ordering> {
        impl_point_cmp_line!(self, other)
    }
}

impl<N: RealField> InstinctOrd<(Point3<N>, Point3<N>)> for Point3<N> {
    fn instinct_cmp(&self, other: &(Point3<N>, Point3<N>)) -> Option<Ordering> {
        impl_point_cmp_line!(self, other)
    }
}

impl<N: RealField> InstinctOrd<(&Point3<N>, &Vector3<N>)> for Point3<N> {
    fn instinct_cmp(&self, other: &(&Point3<N>, &Vector3<N>)) -> Option<Ordering> {
        let end = other.0 + other.1;
        self.instinct_cmp(&(other.0, &end))
    }
}

impl<N: RealField> InstinctOrd<(Point3<N>, Vector3<N>)> for Point3<N> {
    fn instinct_cmp(&self, other: &(Point3<N>, Vector3<N>)) -> Option<Ordering> {
        let end = other.0 + other.1;
        self.instinct_cmp(&(&other.0, &end))
    }
}

//--------------------------------
/// Compare point and a plane linked by three points
macro_rules! impl_point_cmp_plane {
    ($self: ident, $other: ident) => {
        {
            let v = $self.z;
            let min = $other.0.z.min($other.1.z).min($other.2.z);
            let max = $other.0.z.max($other.1.z).max($other.2.z);
            if v >= max {
                return v.partial_cmp(&max);
            }
            if v <= min {
                return v.partial_cmp(&min);
            }

            let line1 = Unit::new_unchecked($other.1.coords - $other.0.coords);
            let line2 = Unit::new_unchecked($other.2.coords - $other.0.coords);
            let target = Unit::new_unchecked($self.coords - $other.0.coords);

            let norm = line1.cross(&line2);
            let delta = norm.dot(&target);
            if norm.z.is_positive() {
                if delta.is_positive() {
                    Some(Ordering::Greater)
                } else if delta.is_negative() {
                    Some(Ordering::Less)
                } else {
                    Some(Ordering::Equal)
                }
            } else if norm.z.is_negative() {
                if delta.is_positive() {
                    Some(Ordering::Less)
                } else if delta.is_negative() {
                    Some(Ordering::Greater)
                } else {
                    Some(Ordering::Equal)
                }
            } else {
                Some(Ordering::Equal)
            }
        }
    }
}
 
impl<N: RealField> InstinctOrd<(Point3<N>, Point3<N>, Point3<N>)> for Point3<N> {
    fn instinct_cmp(&self, other: &(Point3<N>, Point3<N>, Point3<N>)) -> Option<Ordering> {
        impl_point_cmp_plane!(self, other)
    }
}

impl<N: RealField> InstinctOrd<(&Point3<N>, &Point3<N>, &Point3<N>)> for Point3<N> {
    fn instinct_cmp(&self, other: &(&Point3<N>, &Point3<N>, &Point3<N>)) -> Option<Ordering> {
        impl_point_cmp_plane!(self, other)
    }
}
 
macro_rules! reverse_cmp {
    ($self: ident, $other: ident) => {
        match $other.instinct_cmp($self) {
            Some(res) => {
                match res {
                    Ordering::Greater => Some(Ordering::Less),
                    Ordering::Less => Some(Ordering::Greater),
                    Ordering::Equal => Some(Ordering::Equal),
                }
            }
            None => None
        }
    }
}

//--------------------------------
/// Compare line linked by two points and a point
impl<N: RealField> InstinctOrd<Point3<N>> for (Point3<N>, Point3<N>) {
    fn instinct_cmp(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp!(self, other)
    }
}

impl<N: RealField> InstinctOrd<Point3<N>> for (&Point3<N>, &Point3<N>) {
    fn instinct_cmp(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp!(self, other)
    }
}

impl<N: RealField> InstinctOrd<Point3<N>> for (Point3<N>, Vector3<N>) {
    fn instinct_cmp(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp!(self, other)
    }
}

impl<N: RealField> InstinctOrd<Point3<N>> for (&Point3<N>, &Vector3<N>) {
    fn instinct_cmp(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp!(self, other)
    }
}

//--------------------------------
/// Compare two lines linked by two points for each
macro_rules! impl_line_cmp_line {
    ($self: ident, $other: ident) => {
        {
            let v_min = $self.0.z.min($self.1.z);
            let v_max = $self.0.z.max($self.1.z);
            let min = $other.0.z.min($other.1.z);
            let max = $other.0.z.max($other.1.z);
            if v_min >= max {
                return v_min.partial_cmp(&max);
            }
            if v_max <= min {
                return v_max.partial_cmp(&min);
            }
            let c_0 = $self.0.instinct_cmp($other);
            let c_1 = $self.1.instinct_cmp($other);
            if c_0 == c_1 {
                c_0
            } else {
                Some(Ordering::Equal)
            }
        }
    }
}

macro_rules! repeat_impl_line_cmp_line {
    ($a: ty, $b: ty) => {
        impl<N: RealField> InstinctOrd<($a, $a)> for ($b, $b) {
            fn instinct_cmp(&self, other: &($a, $a)) -> Option<Ordering> {
                impl_line_cmp_line!(self, other)
            }
        }
    };
}

repeat_impl_line_cmp_line!(Point3<N>, Point3<N>);
repeat_impl_line_cmp_line!(&Point3<N>, &Point3<N>);
repeat_impl_line_cmp_line!(&Point3<N>, Point3<N>);
repeat_impl_line_cmp_line!(Point3<N>, &Point3<N>);

//--------------------------------
/// Compare line linked by two points and a point
impl<N: RealField> InstinctOrd<Point3<N>> for (Point3<N>, Point3<N>, Point3<N>) {
    fn instinct_cmp(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp!(self, other)
    }
}

impl<N: RealField> InstinctOrd<Point3<N>> for (&Point3<N>, &Point3<N>, &Point3<N>) {
    fn instinct_cmp(&self, other: &Point3<N>) -> Option<Ordering> {
        reverse_cmp!(self, other)
    }
}


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
        assert_eq!(line.instinct_cmp(&Point3::new(3., 5., 5.99)), Some(Ordering::Greater));
        assert_eq!(line.instinct_cmp(&Point3::new(3., 5., 6.001)), Some(Ordering::Less));
    }

    #[test]
    fn test_cmp_point_plane() {
        let plane = (Point3::new(1., 2., 3.), Point3::new(1.4, 4., 6.), Point3::new(0., 0., 0.));
        let mut point = Point3::new(0., 0., 0.);
        assert_eq!(plane.instinct_cmp(&point), Some(Ordering::Equal));
        point = Point3::new(0.8, 2., 3.);
        assert_eq!(plane.instinct_cmp(&point), Some(Ordering::Equal));
        point = Point3::new(0.8, 2., 2.99);
        assert_eq!(plane.instinct_cmp(&point), Some(Ordering::Greater));
        point = Point3::new(0.8, 2., 3.001);
        assert_eq!(plane.instinct_cmp(&point), Some(Ordering::Less));
    }
 
}



