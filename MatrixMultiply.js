
function MatrixVectorMultiply(matrix , vector){
    const[x,y,z] = vector;
    const w=1;
    const result = [
        matrix[0][0]*x + matrix[0][1]*y + matrix[0][2]*z + matrix[0][3]*w,
        matrix[1][0]*x + matrix[1][1]*y + matrix[1][2]*z + matrix[1][3]*w,        
        matrix[2][0]*x + matrix[2][1]*y + matrix[2][2]*z + matrix[2][3]*w,
        matrix[3][0]*x + matrix[3][1]*y + matrix[3][2]*z + matrix[3][3]*w,
     ]

    return result.map(v => v / result[3]); // Perform perspective divide (normalize by wresut[3])
}


function MatrixMultiply(m1,m2){
    const result = [];
    for(let i=0 ; i<4 ; i++){
        result[i] = [];
        for(let j=0; j<4 ; j++){
            result[i][j] = 0;
            for(let k=0; k<4; k++){
                result[i][j] += m1[i][k] * m2[k][j];
            }

        }
    }
}

const matrix = [
    [1, 0, 0, 5],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
];


const vector = [2, 3, 4];

        const result = MatrixMultiply(matrix, vector);
        console.log(result);