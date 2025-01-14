const canvas = document.getElementById("renderer");
        const ctx = canvas.getContext("2d");

        const width = canvas.width;
        const height = canvas.height;

        // Framebuffer
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

        // Matrices and transformations
        function matrixVectorMultiply(matrix, vector) {
            const [x, y, z] = vector;
            const w = 1;
            return [
                matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3] * w,
                matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3] * w,
                matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3] * w,
                matrix[3][0] * x + matrix[3][1] * y + matrix[3][2] * z + matrix[3][3] * w,
            ];
        }

        function matrixMultiply(m1, m2) {
            const result = [];
            for (let i = 0; i < 4; i++) {
                result[i] = [];
                for (let j = 0; j < 4; j++) {
                    result[i][j] = 0;
                    for (let k = 0; k < 4; k++) {
                        result[i][j] += m1[i][k] * m2[k][j];
                    }
                }
            }
            return result;
        }

        function identityMatrix() {
            return [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
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

        function scalingMatrix(sx, sy, sz) {
            return [
                [sx, 0, 0, 0],
                [0, sy, 0, 0],
                [0, 0, sz, 0],
                [0, 0, 0, 1],
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


/*
        const vertexbuffer = [
            [-1, -1, 0],  // Bottom-left corner
            [1, -1, 0],   // Bottom-right corner
            [1, 1, 0],    // Top-right corner
            [-1, 1, 0]    // Top-left corner
        ];
        
        
        const indexbuffer = [
            [0, 1, 2],  // First triangle (bottom-left, bottom-right, top-right)
            [0, 2, 3]   // Second triangle (bottom-left, top-right, top-left)
        ];
        
        
        
        const colorbuffer = [
            [255, 0, 255, 255],  // Red for vertex 0 (bottom-left)
            [255,255 , 255, 255],  // Green for vertex 1 (bottom-right)
            [255, 0, 0, 255],  // Blue for vertex 2 (top-right)
            [255, 255, 0, 255] // Yellow for vertex 3 (top-left)
        ];

        */
     /*
        // Rendering
        const vertexBuffer = [
            [-0.5, -0.5, 0],
            [0.5, -0.5, 0],
            [0, 0.5, 0],
        ];

        const colorBuffer = [
            [255, 0, 0, 255],
            [0, 255, 0, 255],
            [0, 0, 255, 255],
        ];

        const indexbuffer = [
            [0, 1, 2],  // Triangle formed by vertices 0, 1, and 2
        ];
*/


/*
const vertexBuffer = [
    [-1, -1, -1],  // 0
    [1, -1, -1],   // 1
    [1, 1, -1],    // 2
    [-1, 1, -1],   // 3
    [-1, -1, 1],   // 4
    [1, -1, 1],    // 5
    [1, 1, 1],     // 6
    [-1, 1, 1],    // 7
];

const edges = [
    [0, 1], [1, 2], [2, 3], [3, 0],   // Bottom face
    [4, 5], [5, 6], [6, 7], [7, 4],   // Top face
    [0, 4], [1, 5], [2, 6], [3, 7],   // Connecting edges
];

*/

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



        function drawTriangle(v0, v1, v2, c0, c1, c2) {
            // Sort vertices by Y-coordinate (v0.y <= v1.y <= v2.y)
            if (v1[1] < v0[1]) [v0, v1] = [v1, v0];
            if (v2[1] < v1[1]) [v1, v2] = [v2, v1];
            if (v1[1] < v0[1]) [v0, v1] = [v1, v0];
        
            // Helper function to interpolate X between two points
            function interpolateY(x0, y0, x1, y1, y) {
                return x0 + (x1 - x0) * ((y - y0) / (y1 - y0));
            }
        
            // Helper function to interpolate color between two colors
            function interpolateColor(c0, c1, t) {
                return [
                    Math.round(c0[0] + t * (c1[0] - c0[0])),
                    Math.round(c0[1] + t * (c1[1] - c0[1])),
                    Math.round(c0[2] + t * (c1[2] - c0[2])),
                    Math.round(c0[3] + t * (c1[3] - c0[3])) // Interpolating alpha
                ];
            }
        
            // Rasterize the triangle
            for (let y = Math.ceil(v0[1]); y <= Math.floor(v2[1]); y++) {
                let xLeft, xRight;
                let cLeft, cRight;
        
                if (y < v1[1]) {
                    // Top half of the triangle
                    xLeft = interpolateY(v0[0], v0[1], v1[0], v1[1], y);
                    xRight = interpolateY(v0[0], v0[1], v2[0], v2[1], y);
        
                    cLeft = interpolateColor(c0, c1, (y - v0[1]) / (v1[1] - v0[1]));
                    cRight = interpolateColor(c0, c2, (y - v0[1]) / (v2[1] - v0[1]));
                } else {
                    // Bottom half of the triangle
                    xLeft = interpolateY(v1[0], v1[1], v2[0], v2[1], y);
                    xRight = interpolateY(v0[0], v0[1], v2[0], v2[1], y);
        
                    cLeft = interpolateColor(c1, c2, (y - v1[1]) / (v2[1] - v1[1]));
                    cRight = interpolateColor(c0, c2, (y - v0[1]) / (v2[1] - v0[1]));
                }
        
                // Interpolate colors for each pixel on the scanline
                for (let x = Math.ceil(Math.min(xLeft, xRight)); x <= Math.floor(Math.max(xLeft, xRight)); x++) {
                    // Interpolate color along the horizontal scanline (between xLeft and xRight)
                    const t = (x - Math.min(xLeft, xRight)) / (Math.max(xLeft, xRight) - Math.min(xLeft, xRight));
                    const pixelColor = interpolateColor(cLeft, cRight, t);
                    setPixel(x, y, pixelColor[0], pixelColor[1], pixelColor[2], pixelColor[3]);
                }
            }
        }
        


function drawTriangles(vertexBuffer, indexBuffer, color) {
    for (const [i0, i1, i2] of indexBuffer) {
        const v0 = vertexBuffer[i0];
        const v1 = vertexBuffer[i1];
        const v2 = vertexBuffer[i2];

        // Rasterize the triangle
        drawTriangle(
            [v0[0], v0[1]],
            [v1[0], v1[1]],
            [v2[0], v2[1]],
            color[0], color[1], color[2], color[3]
        );
    }
}



/*
function render() {
    
    clearFramebuffer(0, 0, 0, 255);

    const aspectRatio = width / height;
    const projection = perspectiveProjection(90, aspectRatio, 0.1, 100);

    const translation = translationMatrix(0, 0, -5); // Move the cube back
    const rotation = rotationMatrixY(angle);
    const scaling = scalingMatrix(1.5, 1.5, 1.5); // Adjust scale

    // Combine transformations: scaling → rotation → translation → projection
    const transform = matrixMultiply(projection, matrixMultiply(translation, matrixMultiply(rotation, scaling)));

    // Transform vertices
    const transformedVertices = vertexBuffer.map(vertex =>
        matrixVectorMultiply(transform, [...vertex, 1])
    ).map(v => [
        Math.floor((v[0] / v[3] + 1) * width / 2),
        Math.floor((1 - v[1] / v[3]) * height / 2),
    ]);

    // Draw the cube wireframe by connecting the edges
    edges.forEach(edge => {
        const [startIdx, endIdx] = edge;
        const startVertex = transformedVertices[startIdx];
        const endVertex = transformedVertices[endIdx];
        
        drawLine(startVertex[0], startVertex[1], endVertex[0], endVertex[1], [255, 255, 255, 255]); // White edges
    });

    ctx.putImageData(framebuffer, 0, 0);

    angle += 1; // Rotate the object
    requestAnimationFrame(render);
}


let angle = 0;
render();
*/

/*
function render(){
    ClearCanvas();
    clearFramebuffer(0, 0, 0);  // Optional: Clear framebuffer

    

    // Set the transformation matrices for the world, view, and projection
    const worldMatrix = worldMatrix(0, 0, 0, 0, 0, 0, 1, 1, 1);  // Example: Identity transformation
    const viewMatrix = viewMatrix([0, 0, -5], [0, 0, 0], [0, 1, 0]);  // Camera setup
    const projectionMatrix = perspectiveProjection(90, width / height, 0.1, 1000);

    // Transform and project the vertices
    const transformedVertices = transformAndProject(vertexBuffer, worldMatrix, viewMatrix, projectionMatrix);

    // Loop through the transformed vertices and draw them (you can use a triangle drawing function here)
    for (let i = 0; i < indexBuffer.length; i++) {
        const indices = indexBuffer[i];
        const v0 = transformedVertices[indices[0]];
        const v1 = transformedVertices[indices[1]];
        const v2 = transformedVertices[indices[2]];

        // Here you can draw the triangle with interpolated colors
        drawTriangle(v0, v1, v2, colorBuffer[0], colorBuffer[1], colorBuffer[2]);
    }

    // Update the canvas framebuffer
    ctx.putImageData(framebuffer, 0, 0);
    requestAnimationFrame(render)
}

render();
*/



// main rende rworking
/*
function render() {
    clearFramebuffer(0, 0, 0, 255);

    const aspectRatio = width / height;
    const projection = perspectiveProjection(90, aspectRatio, 0.1, 100);

    const translation = translationMatrix(0, 0, -2);
    const rotation = rotationMatriZ(angle);
    const transform = matrixMultiply(projection, matrixMultiply(translation, rotation));

    const transformedVertices = vertexbuffer.map(vertex =>
        matrixVectorMultiply(transform, [...vertex, 1])
    ).map(v => [
        Math.floor((v[0] / v[3] + 1) * width / 2),
        Math.floor((1 - v[1] / v[3]) * height / 2),
    ]);

    // Pass the individual colors for each vertex
    drawTriangle(
        transformedVertices[0], transformedVertices[1], transformedVertices[2],
        colorbuffer[0], colorbuffer[1], colorbuffer[2]
    );
    


    ctx.putImageData(framebuffer, 0, 0);

    angle += 1;
    requestAnimationFrame(render);
}

let angle = 0;
render();
*/




function render() {
    clearFramebuffer(0, 0, 0, 255);

    const aspectRatio = width / height;
    const projection = perspectiveProjection(90, aspectRatio, 0.1, 100);

    const translation = translationMatrix(0, 0, 3);
    const rotation = rotationMatrixY(angle);
    const scaling = scalingMatrix(0.5, 0.5, 0.5); // Example scaling factors for x, y, z

    // Combine transformations: scaling → rotation → translation → projection
    const transform = matrixMultiply(projection, matrixMultiply(translation, matrixMultiply(rotation, scaling)));

    const transformedVertices = vertexbuffer.map(vertex =>
        matrixVectorMultiply(transform, [...vertex, 1])
    ).map(v => [
        Math.floor((v[0] / v[3] + 1) * width / 2),
        Math.floor((1 - v[1] / v[3]) * height / 2),
    ]);

    // Pass the individual colors for each vertex
    drawTriangle(
        transformedVertices[0], transformedVertices[1], transformedVertices[2],
        colorbuffer[0], colorbuffer[1], colorbuffer[2]
    );
    drawTriangle(
        transformedVertices[0], transformedVertices[2], transformedVertices[3],
        colorbuffer[0], colorbuffer[2], colorbuffer[3]
    );

    ctx.putImageData(framebuffer, 0, 0);

    angle += 1; // Rotate the object
    requestAnimationFrame(render);
}

let angle = 0;
render();



/*


function render() {
    
    ClearCanvas();
    clearFrameBuffer(0, 0, 0);

    

    //const worldMatrix = WorldMatrix(0, 0, 0, 0, 0, 0, 1, 1, 1);
   // rotationangle += speed;
    //speed += 0.1;
    const rotationmatrix = RotationMatrixZ(rotationangle);
    const viewMatrix = ViewMatrix([0, 0, -5], [0, 0, 0], [0, 1, 0]);
    const projectionMatrix = PerspectiveProjection(90, width / height, 0.1, 100);

    
    const transformedVertices = transformAndProject(vertexbuffer, rotationmatrix, viewMatrix, projectionMatrix);

    // Draw all triangles/quads in the scene
    drawanything(transformedVertices, indexbuffer, colorbuffer);

    
    ctx.putImageData(framebuffer, 0, 0);

    rotationangle += 1;
    requestAnimationFrame(render);
}



render();

function render() {
    
    ClearCanvas();
    clearFrameBuffer(0, 0, 0);

    

    //const worldMatrix = WorldMatrix(0, 0, 0, 0, 0, 0, 1, 1, 1);
   // rotationangle += speed;
    //speed += 0.1;
    const rotationmatrix = RotationMatrixZ(rotationangle);
    const viewMatrix = ViewMatrix([0, 0, -5], [0, 0, 0], [0, 1, 0]);
    const projectionMatrix = PerspectiveProjection(90, width / height, 0.1, 100);

    
    const transformedVertices = transformAndProject(vertexbuffer, rotationmatrix, viewMatrix, projectionMatrix);

    // Draw all triangles/quads in the scene
    drawanything(transformedVertices, indexbuffer, colorbuffer);

    
    ctx.putImageData(framebuffer, 0, 0);

    rotationangle += 1;
    requestAnimationFrame(render);
}



render();
*/