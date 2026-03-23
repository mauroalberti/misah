use nalgebra::Point3;

fn main() {
    // Definiamo due punti nello spazio 3D
    let p1 = Point3::new(1.0, 2.0, 3.0);
    let p2 = Point3::new(4.0, 6.0, 3.0);

    // Calcolo della distanza euclidea
    let distance = (p2 - p1).norm();

    println!("Distance = {}", distance);
}
