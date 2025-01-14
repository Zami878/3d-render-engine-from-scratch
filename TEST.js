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



function RotationMatrixZ(angle){
    const rad = (angle* Math.PI)/180;
    return[
        [Math.cos(rad),-Math.sin(rad),0,0],
        [Math.sin(rad), Math.cos(rad), 0, 0],
        [0,0,1,0],
        [0, 0, 0, 1],
    ]
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


const faces = [
    [0, 1, 2, 3], 
    [4, 5, 6, 7], 
    [0, 1, 5, 4], 
    [1, 2, 6, 5], 
    [2, 3, 7, 6], 
    [3, 0, 4, 7], 
];


const faceColors = [
    [255, 0, 0],   
    [0, 255, 0],   
    [0, 0, 255],   
    [255, 255, 0], 
    [255, 0, 255], 
    [0, 255, 255], 
];

function drawLineH(y, xStart, xEnd, r, g, b) {
    for (let x = Math.round(xStart); x <= Math.round(xEnd); x++) {
        setPixel(x, y, r, g, b);
    }
}

function fillTriangle(v0, v1, v2, r, g, b) {
    const vertices = [v0, v1, v2].sort((a, b) => a[1] - b[1]);
    const [p0, p1, p2] = vertices;

    const edge1Slope = (p1[0] - p0[0]) / (p1[1] - p0[1] || 1);
    const edge2Slope = (p2[0] - p0[0]) / (p2[1] - p0[1] || 1);
    const edge3Slope = (p2[0] - p1[0]) / (p2[1] - p1[1] || 1);

    let xStart = p0[0];
    let xEnd = p0[0];

    for (let y = Math.round(p0[1]); y <= Math.round(p1[1]); y++) {
        drawHorizontalLine(y, xStart, xEnd, r, g, b);
        xStart += edge1Slope;
        xEnd += edge2Slope;
    }

    xStart = p1[0];

    for (let y = Math.round(p1[1]); y <= Math.round(p2[1]); y++) {
        drawHorizontalLine(y, xStart, xEnd, r, g, b);
        xStart += edge3Slope;
        xEnd += edge2Slope;
    }
}

// Render the cube
let angle = 30;

function render() {
    clearFramebuffer(0, 0, 0, 255);

    const aspectRatio = width / height;
    const projection = perspectiveProjection(90, aspectRatio, 0.1, 1000);

    const translation = translationMatrix(0, 0, 10);
    const rotation = RotationMatrixZ(angle);
    const transform = matrixMultiply(matrixMultiply(projection, translation), rotation);

    const transformedVertices = vertexBuffer.map(vertex => {
        const v = matrixVectorMultiply(transform, [...vertex, 1]);
        return [
            Math.floor((v[0] / v[3]) * width / 2 + width / 2),
            Math.floor((-v[1] / v[3]) * height / 2 + height / 2),
            v[2] / v[3],
        ];
    });

    faces.forEach((face, faceIdx) => {
        const [v0, v1, v2, v3] = face.map(idx => transformedVertices[idx]);
        const [r, g, b] = faceColors[faceIdx];

        fillTriangle(v0, v1, v2, r, g, b);
        fillTriangle(v2, v3, v0, r, g, b);
    });

    ctx.putImageData(framebuffer, 0, 0);

    angle += 1;
    requestAnimationFrame(render);
}

render();
