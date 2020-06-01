use rand::prelude::*;
use rand::distributions::Standard;
use std::iter;
use itertools::multizip;

pub struct Matrix {
    matrix: Vec<f64>,
    size: usize,
}

impl Matrix {
    pub fn new(size: usize) -> Matrix {
        Matrix {
            matrix: iter::repeat(0.0).take(size * size).collect(),
            size: size,
        }
        
    }
    pub fn random(size: usize) -> Matrix {
        let rng = rand::thread_rng();
        let matrix = rng.sample_iter(Standard).take(size * size).collect();
        Matrix {
            matrix: matrix,
            size: size,
        }
        
    }
    pub fn from_vec(elements: Vec<f64>, size: usize) -> Matrix {
        Matrix {
            matrix: elements,
            size: size,
        }
    }
    pub fn pad(&self) -> Matrix {
        let new_size = self.size.next_power_of_two();
        let mut res = Matrix::new(new_size);
        for i in 0..self.size {
            for j in 0..self.size {
                res.set(i, j, self.get(i, j));
            }
        }
        res
    }
    pub fn print(&self) {
        for i in 0..self.size {
            for j in 0..self.size {
                print!("{:.3} ", self.get(i, j));
            }
            println!();
        }
    }
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.matrix[i * self.size + j]
    }
    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        self.matrix[i * self.size + j] = val;
    }
    pub fn get_submatrices(&self) -> [Matrix; 4] {
        let subsize = self.size / 2;
        let mut a = Vec::with_capacity(subsize * subsize);
        let mut b = Vec::with_capacity(subsize * subsize);
        let mut c = Vec::with_capacity(subsize * subsize);
        let mut d = Vec::with_capacity(subsize * subsize);
        for i in 0..subsize {
            for j in 0..subsize {
                a.push(self.get(i, j));
                b.push(self.get(i, j + subsize));
                c.push(self.get(i + subsize, j));
                d.push(self.get(i + subsize, j + subsize));
            }
        }

        let a = Matrix::from_vec(a, subsize);
        let b = Matrix::from_vec(b, subsize);
        let c = Matrix::from_vec(c, subsize);
        let d = Matrix::from_vec(d, subsize);
        [a, b, c, d]
    }
    pub fn from_submatrices(
        c1: &mut Matrix,
        c2: &mut Matrix,
        c3: &mut Matrix,
        c4: &mut Matrix
    ) -> Matrix {
        let subsize = c1.size;
        let size = subsize * 2;
        let mut res = Vec::with_capacity(size * size);
        for _ in 0..subsize {
            res.append(&mut c1.matrix.drain(0..subsize).collect());
            res.append(&mut c2.matrix.drain(0..subsize).collect());
        }
        for _ in 0..subsize {
            res.append(&mut c3.matrix.drain(0..subsize).collect());
            res.append(&mut c4.matrix.drain(0..subsize).collect());
        }

        Matrix::from_vec(res, subsize * 2)
    }
    pub fn get_trimmed(&self, new_size: usize) -> Matrix {
        let mut res = Vec::with_capacity(new_size * new_size);
        for i in 0..new_size {
            for j in 0..new_size {
                res.push(self.get(i, j));
            }
        }
        Matrix::from_vec(res, new_size)
    }
}

pub fn add(m1: &Matrix, m2: &Matrix) -> Matrix {
    let mut res = Vec::with_capacity(m1.size * m1.size);
    for (a, b) in m1.matrix.iter().zip(m2.matrix.iter()) {
        res.push(a + b);
    }
    Matrix::from_vec(res, m1.size)
}

pub fn subtract(m1: &Matrix, m2: &Matrix) -> Matrix {
    let mut res = Vec::with_capacity(m1.size * m1.size);
    for (a, b) in m1.matrix.iter().zip(m2.matrix.iter()) {
        res.push(a - b);
    }
    Matrix::from_vec(res, m1.size)
}

pub fn multiply(m1: &Matrix, m2: &Matrix) -> Matrix {
    if m1.size != m2.size {
        panic!();
    }
    let size = m1.size;
    let mut result = Vec::with_capacity(size * size);

    for i in 0..size {
        for j in 0..size {
            let mut sum = 0.0;
            for k in 0..size {
                sum += m1.get(i, k) * m2.get(k, j);
            }
            result.push(sum);
        }
    }
    Matrix::from_vec(result, size)
}
pub fn strassen(matrix_a: &Matrix, matrix_b: &Matrix) -> Matrix {
    if matrix_a.size != matrix_b.size {
        panic!();
    }
    let size = matrix_a.size;
    let m1 = matrix_a.pad();
    let m2 = matrix_b.pad();

    let res_untrimmed = strassen_mult(&m1, &m2);
    let res = res_untrimmed.get_trimmed(size);
    res
}

fn strassen_mult(a: &Matrix, b: &Matrix) -> Matrix {
    if a.size == 2 {
        let (a1, a2, a3, a4) = (a.get(0, 0), a.get(0, 1), a.get(1, 0), a.get(1, 1));
        let (b1, b2, b3, b4) = (b.get(0, 0), b.get(0, 1), b.get(1, 0), b.get(1, 1));
        let m1 = (a1 + a4) * (b1 + b4);
        let m2 = (a3 + a4) * b1;
        let m3 = a1 * (b2 - b4);
        let m4 = a4 * (b3 - b1);
        let m5 = (a1 + a2) * b4;
        let m6 = (a3 - a1) * (b1 + b2);
        let m7 = (a2 - a4) * (b3 + b4);

        let mut result = Matrix::new(2);
        result.set(0, 0, m1 + m4 - m5 + m7);
        result.set(0, 1, m3 + m5);
        result.set(1, 0, m2 + m4);
        result.set(1, 1, m1 - m2 + m3 + m6);
        return result;
    }
    let [a1, a2, a3, a4] = a.get_submatrices();
    let [b1, b2, b3, b4] = b.get_submatrices();

    let m1 = strassen_mult(&add(&a1, &a4), &add(&b1, &b4));
    let m2 = strassen_mult(&add(&a3, &a4), &b1);    
    let m3 = strassen_mult(&a1, &subtract(&b2, &b4));
    let m4 = strassen_mult(&a4, &subtract(&b3, &b1));
    let m5 = strassen_mult(&add(&a1, &a2), &b4);
    let m6 = strassen_mult(&subtract(&a3, &a1), &add(&b1, &b2));
    let m7 = strassen_mult(&subtract(&a2, &a4), &add(&b3, &b4));

    let subsize = a.size / 2;
    let mut c1 = Matrix::new(subsize);
    let mut c2 = Matrix::new(subsize);
    let mut c3 = Matrix::new(subsize);
    let mut c4 = Matrix::new(subsize);

    // c1
    for (c, e1, e4, e5, e7) in multizip((
        c1.matrix.iter_mut(),
        m1.matrix.iter(),
        m4.matrix.iter(),
        m5.matrix.iter(),
        m7.matrix.iter()
    )) {
        *c = e1 + e4 - e5 + e7;
    }
    // c2
    for (c, e3, e5) in multizip((
        c2.matrix.iter_mut(),
        m3.matrix.iter(),
        m5.matrix.iter()
    )) {
        *c = e3 + e5;
    }
    // c3
    for (c, e2, e4) in multizip((
        c3.matrix.iter_mut(),
        m2.matrix.iter(),
        m4.matrix.iter()
    )) {
        *c = e2 + e4;
    }
    // c4
    for (c, e1, e2, e3, e6) in multizip((
        c4.matrix.iter_mut(),
        m1.matrix.iter(),
        m2.matrix.iter(),
        m3.matrix.iter(),
        m6.matrix.iter()
    )) {
        *c = e1 - e2 + e3 + e6;
    }
    Matrix::from_submatrices(&mut c1, &mut c2, &mut c3, &mut c4)
}

pub fn check(matrix_a: &Matrix, matrix_b: &Matrix) -> bool {
    if matrix_a.size != matrix_b.size {
        println!("Matrix size doesn't match");
        return false;
    }
    for (a, b) in matrix_a.matrix.iter().zip(matrix_b.matrix.iter()) {
        if (a * 1000.0).round() != (b * 1000.0).round() {
            println!("Elements don't match: {} and {}", a, b);
            return false;
        }
    }
    true
}
