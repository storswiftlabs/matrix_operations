use dep::std::println;
mod quantize;

struct Matrix<N> {
    values: [Field; N],
    shape: [Field; 2],
    zero_point: Field,
    scale_reciprocal: Field,
}

impl<N> Matrix<N> {
    fn new(values: [Field; N], shape: [Field; 2]) -> Self {
        Self { values, shape, zero_point: 0, scale_reciprocal: 1}
    }

    fn new_quantized(values: [Field; N], shape: [Field; 2], zero_point: Field, scale_reciprocal: Field) -> Self {
        Self { values, shape, zero_point, scale_reciprocal}
    }

    fn add(self: Matrix<N>, matrix: Self) -> Self {
        if (self.zero_point==0) & (self.scale_reciprocal==1) & (matrix.zero_point==0) & (matrix.scale_reciprocal==1) {
            let values = add(self.values, matrix.values);
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
        else {
            let mut values = [0; N];
            let matrix1_values = self.values;
            let matrix2_values = matrix.values;
            for i in 0..matrix1_values.len(){
                values[i] = quantize::add(matrix1_values[i], self.zero_point, self.scale_reciprocal, matrix2_values[i], matrix.zero_point, matrix.scale_reciprocal, self.zero_point, self.scale_reciprocal);
            }
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
    }

    fn sub(self: Matrix<N>, matrix: Self) -> Self {
        if (self.zero_point==0) & (self.scale_reciprocal==1) & (matrix.zero_point==0) & (matrix.scale_reciprocal==1) {
            let values = sub(self.values, matrix.values);
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
        else {
            let mut values = [0; N];
            let matrix1_values = self.values;
            let matrix2_values = matrix.values;
            for i in 0..matrix1_values.len(){
                values[i] = quantize::sub(matrix1_values[i], self.zero_point, self.scale_reciprocal, matrix2_values[i], matrix.zero_point, matrix.scale_reciprocal, self.zero_point, self.scale_reciprocal);
            }
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
    }

    fn hadamard_product(self, matrix: Self) -> Self {
        if (self.zero_point==0) & (self.scale_reciprocal==1) & (matrix.zero_point==0) & (matrix.scale_reciprocal==1) {
            let values = hadamard_product(self.values, matrix.values);
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
        else {
            let mut values = [0; N];
            let matrix1_values = self.values;
            let matrix2_values = matrix.values;
            for i in 0..matrix1_values.len(){
                values[i] = quantize::mul(matrix1_values[i], self.zero_point, self.scale_reciprocal, matrix2_values[i], matrix.zero_point, matrix.scale_reciprocal, self.zero_point, self.scale_reciprocal);
            }
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
    }

    fn mul(self, matrix: Self, mut result: Self) -> Self{
        if (self.zero_point==0) & (self.scale_reciprocal==1) & (matrix.zero_point==0) & (matrix.scale_reciprocal==1) {
            result.values = mul(self.values, matrix.values, self.shape, matrix.shape, result.values);
            result
        }
        else {
            assert(self.shape[1]==matrix.shape[0]);
            for i in 0..self.shape[0]{
                for j in 0..matrix.shape[1]{
                    let mut element = result.zero_point;
                    for k in 0..self.shape[1]{
                        element = quantize::add(quantize::mul(self.values[i * self.shape[1] + k], self.zero_point, self.scale_reciprocal, matrix.values[k * matrix.shape[1] + j], matrix.zero_point, matrix.scale_reciprocal, result.zero_point, result.scale_reciprocal), result.zero_point, result.scale_reciprocal, element, result.zero_point, result.scale_reciprocal, result.zero_point, result.scale_reciprocal);
                    }
                    result.values[i*matrix.shape[1]+j] = element;
                }
            }
            result
        }
    }

    fn scalar_mul(self, scalar: Field) -> Self {
        if (self.zero_point==0) & (self.scale_reciprocal==1) {
            let values = scalar_mul(self.values, scalar);
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
        else {
            let mut values = [0; N];
            for i in 0..N{
                values[i] = quantize::mul(scalar, self.zero_point, self.scale_reciprocal, self.values[i], self.zero_point, self.scale_reciprocal, self.zero_point, self.scale_reciprocal);
            }
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
    }

    fn scalar_div(self, scalar: Field) -> Self {
        if (self.zero_point==0) & (self.scale_reciprocal==1) {
            let values = scalar_div(self.values, scalar);
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
        else {
            let mut values = [0; N];
            for i in 0..N{
                values[i] = quantize::div(self.values[i], self.zero_point, self.scale_reciprocal, scalar, self.zero_point, self.scale_reciprocal, self.zero_point, self.scale_reciprocal);
            }
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
    }

    fn transpose(self) -> Self {
        let values = transpose(self.values, self.shape);
        Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
    }

    fn determinant(self) -> Field{
        if (self.zero_point==0) & (self.scale_reciprocal==1) {
            let det = determinant(self.values, self.shape);
            det
        }
        else {
            assert(self.shape[0] == self.shape[1]);
            if self.shape[0] == 2{
                self.values[0] * self.values[3] - self.values[1] * self.values[2]
            }
            else if self.values[0] == 1{
                self.values[0]
            }
            else {
                let mut det = self.zero_point;
                let mut sign = 1;
                for i in 0..self.values[0]{
                    let mut minor = [0; N];
                    minor = get_minor(self.values, 0, i, self.shape, minor);
                    let matrix = Self {values: minor, shape: [self.shape[0] - 1, self.shape[1] - 1], zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal};
                    det = quantize::add(det, self.zero_point,self.scale_reciprocal, quantize::mul(quantize::mul(sign, 0, 1, self.values[i], self.zero_point,self.scale_reciprocal, self.zero_point,self.scale_reciprocal), self.zero_point,self.scale_reciprocal, matrix.determinant(), self.zero_point,self.scale_reciprocal, self.zero_point,self.scale_reciprocal), self.zero_point,self.scale_reciprocal, self.zero_point,self.scale_reciprocal);
                    sign *= -1;
                }
                det
            }
        }
    }

    fn adjoint(self) -> Self{
        if (self.zero_point==0) & (self.scale_reciprocal==1) {
            let values = adjoint(self.values, self.shape);
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
        else {
            assert(self.shape[0] == self.shape[1]);
            let mut values = [0; N];
            for i in 0..self.shape[0]{
                for j in 0..self.shape[1]{
                    let mut minor = [0; N];
                    minor = get_minor(self.values, i, j, self.shape, minor);
                    let matrix = Self {values: minor, shape: [self.shape[0] - 1, self.shape[1] - 1], zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal};
                    let mut sign = 1;
                    for _ in 0..i+j {
                        sign *= -1;
                    }
                    let mut cofactor = quantize::mul(sign, 0, 1, matrix.determinant(), self.zero_point,self.scale_reciprocal, self.zero_point,self.scale_reciprocal);
                    values[i*self.shape[1]+j] = cofactor;
                }
            }
            values = transpose(values, self.shape);
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
    }

    fn inverse(self) -> Self{
        if (self.zero_point==0) & (self.scale_reciprocal==1) {
            let values = inverse(self.values, self.shape);
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
        else {
            assert(self.shape[0] == self.shape[1]);
            let det = self.determinant();
            assert(det != 0);
            let adjoint = self.adjoint();
            let mut values = [0; N];
            for i in 0..self.shape[0]{
                for j in 0..self.shape[1]{
                    values[i * self.shape[1] + j] = quantize::div(adjoint.values[i * self.shape[1] + j], self.zero_point,self.scale_reciprocal, det, self.zero_point,self.scale_reciprocal, self.zero_point,self.scale_reciprocal);
                }
            }
            Self {values, shape: self.shape, zero_point: self.zero_point, scale_reciprocal: self.scale_reciprocal}
        }
    }

    fn print(self) {
        dep::std::println(self.values);
    }
}

fn add<M>(matrix1: [Field; M], matrix2: [Field; M]) -> [Field; M]{
    let mut result = [0; M];
    for i in 0..matrix1.len(){
        result[i] = (matrix1[i] + matrix2[i]);
    }
    result
}

fn sub<M>(matrix1: [Field; M], matrix2: [Field; M]) -> [Field; M]{
    let mut result = [0; M];
    for i in 0..matrix1.len(){
        result[i] = (matrix1[i] - matrix2[i]);
    }
    result
}

fn hadamard_product<M>(matrix1: [Field; M], matrix2: [Field; M]) -> [Field; M]{
    let mut result = [0; M];
    for i in 0..matrix1.len(){
        result[i] = (matrix1[i] * matrix2[i]);
    }
    result
}

fn mul<M, N>(matrix1: [Field; M], matrix2: [Field; M], matrix1_shape: [comptime Field; 2], matrix2_shape: [comptime Field; 2], mut result: [Field; N]) -> [Field; N]{
    assert(matrix1_shape[1]==matrix2_shape[0]);
    for i in 0..matrix1_shape[0]{
        for j in 0..matrix2_shape[1]{
            let mut element = 0;
            for k in 0..matrix1_shape[1]{
                element += matrix1[i * matrix1_shape[1] + k] * matrix2[k * matrix2_shape[1] + j];
            }
            result[i*matrix2_shape[1]+j] = element;
        }
    }
    result
}

fn scalar_mul<M>(matrix: [Field; M], scalar: Field) -> [Field; M]{
    let mut result = [0; M];
    for i in 0..matrix.len(){
        result[i] = scalar*matrix[i];
    }
    result
}

fn scalar_div<M>(matrix: [Field; M], scalar: Field) -> [Field; M]{
    let mut result = [0; M];
    for i in 0..matrix.len(){
        result[i] = matrix[i]/scalar;
    }
    result
}

fn transpose<M>(matrix: [Field; M], shape: [comptime Field; 2]) -> [Field; M]{
    let mut result = [0; M];
    for i in 0..shape[0]{
        for j in 0..shape[1]{
            result[j * shape[0] + i] = matrix[i * shape[1] + j];
        }
    }
    result
}

fn get_minor<M, N>(matrix: [Field; M], row: Field, col: Field, shape: [comptime Field; 2], mut minor: [Field; N]) -> [Field; N]{
    let mut r: Field = -1;
    for i in 0..shape[0]{
        if i != row{
            r+=1;
            let mut c: Field = -1;
            for j in 0..shape[1]{
                if j != col{
                    c+=1;
                    minor[r*(shape[1]-1)+c]=(matrix[i * shape[1] + j]);
                }
            }
        }
    }
    minor
}

fn determinant<M>(matrix: [Field; M], shape: [comptime Field; 2]) -> Field{
    assert(shape[0] == shape[1]);
    
    if shape[0] == 2{
        matrix[0] * matrix[3] - matrix[1] * matrix[2]
    }
    else if shape[0] == 1{
        matrix[0]
    }
    else {
        let mut det = 0;
        let mut sign = 1;
        for i in 0..shape[0]{
            let mut minor = [0; M];
            minor = get_minor(matrix, 0, i, shape, minor);
            det += sign * matrix[i] * determinant(minor, [shape[0] - 1, shape[1] - 1]);
            sign *= -1;
        }
        det
    }
}

fn adjoint<M>(matrix: [Field; M], shape: [comptime Field; 2]) -> [Field; M]{
    assert(shape[0] == shape[1]);
    
    let mut adjoint = [0; M];
    for i in 0..shape[0]{
        for j in 0..shape[1]{
            let mut minor = [0; M];
            minor = get_minor(matrix, i, j, shape, minor);
            let mut sign = 1;
            for _ in 0..i+j {
                sign *= -1;
            }
            let mut cofactor = sign * determinant(minor, [shape[0] - 1, shape[1] - 1]);
            adjoint[i*shape[1]+j] = cofactor;
        }
    }
    adjoint = transpose(adjoint, shape);
    adjoint
}

fn inverse<M>(matrix: [Field; M], shape: [comptime Field; 2]) -> [Field; M]{
    assert(shape[0] == shape[1]);
    
    let det = determinant(matrix, shape);
    assert(det != 0);
    
    let adjoint = adjoint(matrix, shape);
    
    let mut inverse = [0; M];
    for i in 0..shape[0]{
        for j in 0..shape[1]{
            inverse[i * shape[1] + j] = adjoint[i * shape[1] + j] / det;
        }
    }
    inverse
}

fn main() {
    let matrix1 = Matrix::new([1; 9], [3, 3]);
    let matrix2 = Matrix::new([2; 9], [3, 3]);
    let result = matrix1.add(matrix2);  
    result.print();
}

#[test]
fn test_add() {
    let matrix1 = [1; 9];
    let matrix2 = [2; 9];
    let result = add(matrix1, matrix2);
    println(result);
}

#[test]
fn test_add1() {
    let matrix1 = Matrix::new([1; 9], [3, 3]);
    let matrix2 = Matrix::new([2; 9], [3, 3]);
    let result = matrix1.add(matrix2);  
    result.print();
}

#[test]
fn test_add2() {
    let matrix1 = Matrix::new_quantized([1; 9], [3, 3], 0, 1);
    let matrix2 = Matrix::new_quantized([2; 9], [3, 3], 0, 2);
    let result = matrix1.add(matrix2);  // {"scale_reciprocal":"0x01","shape":["0x03","0x03"],"values":["0x05","0x05","0x05","0x05","0x05","0x05","0x05","0x05","0x05"],"zero_point":"0x00"}
    println(result);
}

#[test]
fn test_sub() {
    let matrix1 = [3; 9];
    let matrix2 = [2; 9];
    let result = sub(matrix1, matrix2);
    println(result);
}

#[test]
fn test_hadamard_product() {
    let matrix1 = [3; 9];
    let matrix2 = [2; 9];
    let result = hadamard_product(matrix1, matrix2);
    println(result);
}

#[test]
fn test_mul() {
    let matrix1 = [1, 2, 3, 4, 5, 6];
    let matrix2 = [7, 8, 9, 10, 11, 12];
    let matrix1_shape = [2,3];
    let matrix2_shape = [3,2];
    let mut result = [0; 4];
    result = mul(matrix1, matrix2, matrix1_shape, matrix2_shape, result);
    println(result);    //  [58, 64, 139, 154] ['0x3a', '0x40', '0x8b', '0x9a']
}

#[test]
fn test_scalar_mul() {
    let matrix = [1, 2, 3, 4, 5, 6];
    let scalar = 2;
    let result = scalar_mul(matrix, scalar);
    println(result);
}

#[test]
fn test_scalar_div() {
    let matrix = [2, 4, 6];
    let scalar = 2;
    let result = scalar_div(matrix, scalar);
    println(result);
}

#[test]
fn test_transpose() {
    let matrix = [1, 2, 3, 4, 5, 6];
    let shape = [2, 3];
    let result = transpose(matrix, shape);
    println(result);    // [1, 4, 2, 5, 3, 6] ['0x1', '0x4', '0x2', '0x5', '0x3', '0x6']
}

#[test]
fn test_determinant() {
    let matrix = [0, 1, 2, 4, 5, 6, 7, 8, 0];
    let shape = [3, 3];
    let result = determinant(matrix, shape);
    println(result);    // 36 '0x24'
}

#[test]
fn test_adjoint() {
    let matrix = [0, 1, 2, 4, 5, 6, 7, 8, 0];
    let shape = [3, 3];
    let result = adjoint(matrix, shape);
    println(result);    // ['-0x30', '0x10', '-0x4', '0x2a', '-0xe', '0x8', '-0x3', '0x7', '-0x4']
}

#[test]
fn test_inverse() {
    let matrix = [0, 2, 3, 4, 5, 6, 7, 8, 9];
    let shape = [3, 3];
    let result = inverse(matrix, shape);
    println(result);    // ['-0x30', '0x10', '-0x4', '0x2a', '-0xe', '0x8', '-0x3', '0x7', '-0x4']
}
