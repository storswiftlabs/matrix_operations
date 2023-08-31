# Matrix Operations Library

The Matrix Operations Library provides functionality for performing various matrix operations, including addition, subtraction, Hadamard product, matrix multiplication, scalar multiplication, scalar division, transposition, determinant computation, adjoint computation, and matrix inversion. This library is designed to handle both regular and quantized matrices.

## Features

### Matrix Struct

The `Matrix` struct represents a matrix and contains the following attributes:

- `values`: An array holding the matrix values.
- `shape`: An array defining the shape of the matrix.
- `zero_point`: The zero point for quantization.
- `scale`: The scaling factor for quantization.

### Constructors

The `Matrix` struct provides two constructors:

1. `new(values: [Field; N], shape: [comptime Field; 2]) -> Self`: Constructs a new matrix with default quantization parameters.
2. `new_quantized(values: [Field; N], shape: [comptime Field; 2], zero_point: Field, scale: Field) -> Self`: Constructs a new quantized matrix with custom quantization parameters.

### Operations

The following operations are available for matrices:

- Addition: `add(matrix: Self) -> Self`
- Subtraction: `sub(matrix: Self) -> Self`
- Hadamard Product: `hadamard_product(matrix: Self) -> Self`
- Matrix Multiplication: `mul(matrix: Self, result: Self) -> Self`
- Scalar Multiplication: `scalar_mul(scalar: Field) -> Self`
- Scalar Division: `scalar_div(scalar: Field) -> Self`
- Transposition: `transpose() -> Self`
- Determinant Computation: `determinant() -> Field`
- Adjoint Computation: `adjoint() -> Self`
- Matrix Inversion: `inverse() -> Self`

### Printing

The method `print()` allows printing the matrix to the console.

## Usage

1. Import the `Matrix` struct and relevant functions.
2. Create matrices using the provided constructors.
3. Perform desired matrix operations using the available methods.

## Examples

```rust
#[test]
fn test_add() {
    let matrix1 = [1; 9];
    let matrix2 = [2; 9];
    let result = add(matrix1, matrix2);
    println(result);   // ["0x03","0x03","0x03","0x03","0x03","0x03","0x03","0x03","0x03"]
    assert(result==[3; 9]);
}

#[test]
fn test_transpose() {
    let matrix1 = Matrix::new([1, 2, 3, 4, 5, 6], [2, 3]);
    let result = matrix1.transpose(); 
    println(result);     // {"scale":"0x01","shape":["0x03","0x02"],"values":["0x01","0x04","0x02","0x05","0x03","0x06"],"zero_point":"0x00"}
    let transpose = Matrix::new([1, 4, 2, 5, 3, 6], [3, 2]);
    assert(result.shape==transpose.shape);
    assert(result.values==transpose.values);
}

#[test]
fn test_mul() {
    let matrix1 = Matrix::new_quantized([999900, 1000200, 1000300, 1000400, 1000500, 999400], [2, 3], 1000000, 100);    // -quantize-> [-1, 2, 3, 4, 5, -6]
    let matrix2 = Matrix::new_quantized([1000700, 1000800, 1000900, 1001000, 1001100, 998800], [3, 2], 1000000, 100);    // -quantize-> [7, 8, 9, 10, 11, -12]
    let mut result = Matrix::new_quantized([0; 4], [2, 2], 1000000, 100);
    result = matrix1.mul(matrix2, result);  
    result.print();     // ["0x0f5370","0x0f38e0","0x0f44fc","0x0f7e68"] -quantize-> [44, -24, 7, 154]
    assert(result.values==[1004400, 997600, 1000700, 1015400]);
}
```

## Dependencies

This library uses [quantized_arithmetic](https://github.com/storswiftlabs/quantized_arithmetic/tree/main) for quantized arithmetic operations.

## Nargo Version Information

- Nargo Version: 0.9.0
- Git Version Hash: 2232a387177f827d14f7c0b7ac1f3e5bb6d957d9
- Is Dirty: false

## License

This Matrix Operations Library is provided under the [Apache License, Version 2.0](LICENSE).

Feel free to contribute, report issues, and suggest improvements!
