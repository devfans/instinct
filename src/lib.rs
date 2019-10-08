use nalgebra as na;
use nalgebra::*;
use std::cmp::{Ord, PartialOrd, Ordering};


impl<N> PartialOrd for Point3<N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.z.partial_cmp(other.z)
    }
}
