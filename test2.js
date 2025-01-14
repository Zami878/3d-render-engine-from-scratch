const canvas = document.getElementById("renderer");
const ctx = canvas.getContext("2d");

const width = canvas.width;
const height = canvas.height;

const framebuffer = ctx.createImageData(width, height);

function setPixel(x, y, r, g, b, a = 255) {
    if (x < 0 || y < 0 || x >= width || y >= height) return;
    const index = (y * width + x) * 4;
    framebuffer.data[index] = r;
    framebuffer.data[index + 1] = g;
    framebuffer.data[index + 2] = b;
    framebuffer.data[index + 3] = a;
}

function clearFramebuffer(r, g, b, a = 255) {
    for (let i = 0; i < framebuffer.data.length; i += 4) {
        framebuffer.data[i] = r;
        framebuffer.data[i + 1] = g;
        framebuffer.data[i + 2] = b;
        framebuffer.data[i + 3] = a;
    }
}

function matrixVectorMultiply(matrix, vector) {
    const [x, y, z, w = 1] = vector;
    return [
        matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3] * w,
        matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3] * w,
        matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3] * w,
        matrix[3][0] * x + matrix[3][1] * y + matrix[3][2] * z + matrix[3][3] * w,
    ];
}

function perspectiveProjection(fovy, aspectRatio, znear, zfar) {
    const h = 1 / Math.tan((fovy * 0.5) * (Math.PI / 180));
    const w = h / aspectRatio;
    const a = zfar / (zfar - znear);
    const b = (-znear * zfar) / (zfar - znear);

    return [
        [w, 0, 0, 0],
        [0, h, 0, 0],
        [0, 0, a, b],
        [0, 0, -1, 0],
    ];
}

function translationMatrix(tx, ty, tz) {
    return [
        [1, 0, 0, tx],
        [0, 1, 0, ty],
        [0, 0, 1, tz],
        [0, 0, 0, 1],
    ];
}

function rotationMatrixY(angle) {
    const rad = (angle * Math.PI) / 180;
    return [
        [Math.cos(rad), 0, Math.sin(rad), 0],
        [0, 1, 0, 0],
        [-Math.sin(rad), 0, Math.cos(rad), 0],
        [0, 0, 0, 1],
    ];
}

function matrixMultiply(m1, m2) {
    const result = Array.from({ length: 4 }, () => Array(4).fill(0));
    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            for (let k = 0; k < 4; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return result;
}

/*
const vertexBuffer = [
    [-2, -1, -1],  
    [2, -1, -1],   
    [2, 1, -1],    
    [-2, 1, -1],   
    [-2, -1, 1],   
    [2, -1, 1],   
    [2, 1, 1],    
    [-2, 1, 1],    
];

const edges = [
    [0, 1], [1, 2], [2, 3], [3, 0],
    [4, 5], [5, 6], [6, 7], [7, 4],
    [0, 4], [1, 5], [2, 6], [3, 7],
];
*/




const vertexBuffer = [
    [-1, -1, -1],  
    [1, -1, -1],   
    [1, 1, -1],    
    [-1, 1, -1],   
    [-1, -1, 1],   
    [1, -1, 1],    
    [1, 1, 1],     
    [-1, 1, 1],    
];

// kind of like ondex buffer 
const edges = [
    [0, 1], [1, 2], [2, 3], [3, 0],   
    [4, 5], [5, 6], [6, 7], [7, 4],   
    [0, 4], [1, 5], [2, 6], [3, 7],  
];



// bersenham all octants
function drawLine(x0, y0, x1, y1, r, g, b, a = 255) {
    const dx = Math.abs(x1 - x0);
    const dy = Math.abs(y1 - y0);
    const sx = x0 < x1 ? 1 : -1;
    const sy = y0 < y1 ? 1 : -1;
    let err = dx - dy;

    while (true) {
        setPixel(x0, y0, r, g, b, a);
        if (x0 === x1 && y0 === y1) break;
        const e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x0 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y0 += sy;
        }
    }
}

let angle = 0;

function render() {
    clearFramebuffer(0, 0, 0, 255);

    const aspectRatio = width / height;
    const projection = perspectiveProjection(90, aspectRatio, 0.1, 1000);
    
    // back in  z axis
    const translation = translationMatrix(0, 0, 5);
    const rotation = rotationMatrixY(angle);
    const transform = matrixMultiply(matrixMultiply(projection, translation), rotation);

    
    const transformedVertices = vertexBuffer.map(vertex => {
        const v = matrixVectorMultiply(transform, [...vertex, 1]);
        return [
            Math.floor((v[0] / v[3]) * width / 2 + width / 2), // ro map correctly [-width/2 , width/2] => [0,width/2]
            Math.floor((-v[1] / v[3]) * height / 2 + height / 2),
        ];
    });

  
    edges.forEach(([startIdx, endIdx]) => {   // edges pair start or emnd pts
        const startVertex = transformedVertices[startIdx];
        const endVertex = transformedVertices[endIdx];

        drawLine(
            startVertex[0],
            startVertex[1],
            endVertex[0],
            endVertex[1],
            255, 255, 255, 255
        );
    });

    ctx.putImageData(framebuffer, 0, 0);

    
    angle += 1;

    requestAnimationFrame(render);
}   


render();
