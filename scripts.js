const canvas = document.getElementById("renderer");

const ctx = canvas.getContext("2d");


function ClearCanvas(){
    ctx.fillStyle = "black";
    ctx.fillRect(0,0,canvas.width , canvas.height);
}



const width = canvas.width;
const height = canvas.height;

//framebuffer
const framebuffer = ctx.createImageData(width,height);

//a pixel in framebuffer
function setPixel(x,y,r,g,b,a = 255) {
    if(x<0 || y<0 || x>= width || y>= height) return;
    const index = (y*width+x)*4;
    framebuffer.data[index] = r;
    framebuffer.data[index+1] = g;
    framebuffer.data[index+2] = b;
    framebuffer.data[index+3] = a;
}


// 2 buffer strategy 

function clearFrameBuffer(r,g,b,a=255){
    for(let i=0; i<framebuffer.data.length; i+=4){
        framebuffer.data[i] = r; 
        framebuffer.data[i+1] = g;
        framebuffer.data[i+2] = b;
        framebuffer.data[i+3] = a;
     }
}

function MatrixVectorMultiply(matrix, vector) {
    console.log("Matrix:", matrix);
    console.log("Vector before transformation:", vector);

    const [x, y, z, w = 1] = vector;  // Ensure w is 1 by default
    
    // Matrix multiplication
    const result = [
        matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3] * w,
        matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3] * w,
        matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3] * w,
        matrix[3][0] * x + matrix[3][1] * y + matrix[3][2] * z + matrix[3][3] * w,
    ];

    console.log("Vector after transformation:", result);
    return result;
}





function MatrixMultiply(m1, m2) {
    if (!m1 || !m2) {
        throw new Error("One or both matrices are undefined");
    }
    
    if (m1.length !== 4 || m2.length !== 4 || m1[0].length !== 4 || m2[0].length !== 4) {
        throw new Error("Both matrices must be 4x4.");
    }

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


/*
const vertexbuffer = [
    [-0.5, -0.5,  0.5], [ 0.5, -0.5,  0.5], [ 0.5,  0.5,  0.5], [-0.5,  0.5,  0.5], // Front face
    [-0.5, -0.5, -0.5], [ 0.5, -0.5, -0.5], [ 0.5,  0.5, -0.5], [-0.5,  0.5, -0.5]  // Back face
];

// Define edges as pairs of indices
const edgebuffer = [
    [0, 1], [1, 2], [2, 3], [3, 0], // Front face
    [4, 5], [5, 6], [6, 7], [7, 4], // Back face
    [0, 4], [1, 5], [2, 6], [3, 7]  // Connecting edges
];
*/


/*
const vertexbuffer = [
    [-1, -1, 0],  // Vertex 0
    [0.5, -1, 0],   // Vertex 1
    [0, 1, 0],      // Vertex 2
];

// Index buffer: Stores the indices of the vertices that form a triangle
const indexbuffer = [
    [0, 1, 2],  // Triangle formed by vertices 0, 1, and 2
];

// Color buffer: Stores the color of each triangle (RGBA format)
const colorbuffer = [
    [255, 0, 0,255],
    [0,0,0,255],
    [0,0,0,255],
];
*/

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
    [255, 0, 0, 255],  // Red for vertex 0 (bottom-left)
    [255, , 0, 255],  // Green for vertex 1 (bottom-right)
    [255, 0, 0, 255],  // Blue for vertex 2 (top-right)
    [255, 0, 0, 255] // Yellow for vertex 3 (top-left)
];

*/

/*

const vertexbuffer = [
    // First Square (Left)
    [-0.1, -0.1,  0], 
    [ 0.1, -0.1,  0], 
    [ 0.1,  0.1,  0], 
    [-0.1,  0.1,  0], 

    // Second Square (Right, shifted by 0.4 units along X-axis)
    [ 0.6, -0.1,  0], 
    [ 0.8, -0.1,  0], 
    [ 0.8,  0.1,  0], 
    [ 0.6,  0.1,  0]
];

const indexbuffer = [
    // First square (Left)
    [0, 1, 2],
    [0, 2, 3],

    // Second square (Right)
    [4, 5, 6],
    [4, 6, 7]
];


const colorbuffer = [
    [255, 0, 0, 255],   // Red for vertex 0
    [0, 255, 0, 255],   // Green for vertex 1
    [0, 0, 255, 255],   // Blue for vertex 2
    [255, 255, 0, 255], // Yellow for vertex 3
    [255, 0, 255, 255], // Purple for vertex 4
    [0, 255, 255, 255], // Cyan for vertex 5
    [255, 165, 0, 255], // Orange for vertex 6
    [255, 255, 255, 255] // White for vertex 7
];

*/

/*

const vertexbuffer = [
    [0, 1, 0],     // Apex
    [-1, -1, -1],  // Base corner 1
    [1, -1, -1],   // Base corner 2
    [1, -1, 1],    // Base corner 3
    [-1, -1, 1]    // Base corner 4
];

const indexbuffer = [
    [0, 1, 2], // First face
    [0, 2, 3], // Second face
    [0, 3, 4], // Third face
    [0, 4, 1], // Fourth face
    [1, 2, 3],  // Base triangle 1
    [3, 4, 1]   // Base triangle 2
];

// Vertex color buffer
const colorbuffer = [
    [255, 0, 0],    // Apex: Red
    [0, 255, 0],    // Base corner 1: Green
    [0, 0, 255],    // Base corner 2: Blue
    [255, 255, 0],  // Base corner 3: Yellow
    [255, 0, 255]   // Base corner 4: Magenta
];

console.log("Vertex Buffer:", vertexbuffer);
console.log("Index Buffer:", indexbuffer);
console.log("Color Buffer:", colorbuffer);

*/

/*
        // Rendering
const vertexbuffer = [
            [-0.5, -0.5, 0],
            [0.5, -0.5, 0],
            [0, 0.5, 0],
        ];

const colorbuffer = [
            [255, 0, 0, 255],
            [0, 255, 0, 255],
            [0, 0, 255, 255],
        ];

 const indexbuffer = [
            // First square (Left)
            [0, 1, 2],
            [0, 2, 3],
 ];

 */

/*
// Vertex buffer for a cube
// Define the cube vertices
const vertexbuffer = [
    // Front face
    [-0.05, -0.05,  0.05], 
    [ 0.05, -0.05,  0.05], 
    [ 0.05,  0.05,  0.05], 
    [-0.05,  0.05,  0.05], 

    // Back face
    [-0.05, -0.05, -0.05], 
    [ 0.05, -0.05, -0.05], 
    [ 0.05,  0.05, -0.05], 
    [-0.05,  0.05, -0.05], 

    // Left face
    [-0.05, -0.05, -0.05], 
    [-0.05, -0.05,  0.05], 
    [-0.05,  0.05,  0.05], 
    [-0.05,  0.05, -0.05], 

    // Right face
    [ 0.05, -0.05, -0.05], 
    [ 0.05, -0.05,  0.05], 
    [ 0.05,  0.05,  0.05], 
    [ 0.05,  0.05, -0.05], 

    // Top face
    [-0.05,  0.05, -0.05], 
    [-0.05,  0.05,  0.05], 
    [ 0.05,  0.05,  0.05], 
    [ 0.05,  0.05, -0.05], 

    // Bottom face
    [-0.05, -0.05, -0.05], 
    [-0.05, -0.05,  0.05], 
    [ 0.05, -0.05,  0.05], 
    [ 0.05, -0.05, -0.05]
];


// Define the indices to draw the triangles
const indexbuffer = [
    // Front face
    [0, 1, 2],
    [0, 2, 3],

    // Back face
    [4, 5, 6],
    [4, 6, 7],

    // Left face
    [8, 9, 10],
    [8, 10, 11],

    // Right face
    [12, 13, 14],
    [12, 14, 15],

    // Top face
    [16, 17, 18],
    [16, 18, 19],

    // Bottom face
    [20, 21, 22],
    [20, 22, 23]
];

// Define the color buffer (yellow for all vertices)
const colorbuffer = new Array(24).fill([255, 255, 0, 255]);  // Yellow color

*/
/*
const vertexbuffer = [
    [-0.5, -0.5,  0.5],
    [ 0.5, -0.5,  0.5],
    [ 0.5,  0.5,  0.5],
    [-0.5,  0.5,  0.5],
    [-0.5, -0.5, -0.5],
    [ 0.5, -0.5, -0.5],
    [ 0.5,  0.5, -0.5],
    [-0.5,  0.5, -0.5]
];

const indexbuffer = [
    [0, 1], [1, 2], [2, 3], [3, 0],  // Front face
    [4, 5], [5, 6], [6, 7], [7, 4],  // Back face
    [0, 4], [1, 5], [2, 6], [3, 7],  // Connecting front and back faces
];

*/
function identityMatrix() {
    return [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ];
}

function translationMatrix(tx,ty,tz){
    return[
        [1,0,0,tx],
        [0,1,0,ty],
        [0,0,1,tz],
        [0,0,0,1],
    ];
}

function RotationMatrixY(angle_y){
    const rad = (angle_y * Math.PI)/180;
    return[
        [Math.cos(rad),0,Math.sin(rad),0],
        [0, 1, 0, 0],
        [-Math.sin(rad), 0, Math.cos(rad), 0],
        [0, 0, 0, 1],
    ]
}

function RotationMatrixX(angle_x){
    const rad = (angle_x * Math.PI)/180;
    return[
        [1,0,0,0],
        [0, Math.cos(rad), -Math.sin(rad), 0],
        [0,Math.sin(rad), Math.cos(rad), 0],
        [0, 0, 0, 1],
    ]
}



function RotationMatrixZ(angle_z){
    const rad = (angle_z * Math.PI)/180;
    return[
        [Math.cos(rad),-Math.sin(rad),0,0],
        [Math.sin(rad), Math.cos(rad), 0, 0],
        [0,0,1,0],
        [0, 0, 0, 1],
    ]
}

function ROTATION(angle_x , angle_y , angle_z){
    const Rotation_X = RotationMatrixX(angle_x);
    const Rotation_Y = RotationMatrixY(angle_y);
    const Rotation_Z = RotationMatrixZ(angle_z);

        // Check if the individual matrices are 4x4
        if (Rotation_X.length !== 4 || Rotation_Y.length !== 4 || Rotation_Z.length !== 4) {
            console.error("Error: Rotation matrices must be 4x4.");
            return null;  // Or throw an error
        }

    
const Rotation_zy = MatrixMultiply(Rotation_Z , Rotation_Y);
return MatrixMultiply(Rotation_zy, Rotation_X)
}



function ScalingMatrix(sx, sy, sz) {
    return [
        [sx, 0, 0, 0],
        [0, sy, 0, 0],
        [0, 0, sz, 0],
        [0, 0, 0, 1]
    ];
}

function WorldMatrix(tx,ty,tz,angle_x,angle_y,angle_z,sx,sy,sz){

    const Translation = translationMatrix(tx,ty,tz);
    const Rotation = ROTATION(angle_x,angle_y , angle_z);
    const Scaling = ScalingMatrix(sx,sy,sz);

    // world = translATION * rotation * scale
    return MatrixMultiply(MatrixMultiply(Translation ,Rotation),Scaling);

}


function ViewMatrix(eye , center , up){
    //camera local axes
    const[ex,ey,ez] = eye;       //camera posn
    const[cx,cy,cz] = center; // look at  // since camera facing z so correspons to z //
    const[ux,uy,uz] = up  //orientation  , standing upright this is vectoir from head to sky
// right vector correspon x axis 

//forward vector

const xcomponent = cx - ex;
const ycomponent = cy - ey;
const zcomponent = cz - ez;
const distance = Math.sqrt(xcomponent*xcomponent + ycomponent*ycomponent + zcomponent*zcomponent);
const forward = [xcomponent/distance, ycomponent/distance , zcomponent/distance];  //unit verctor dirn important

///MOVE WORLD INSTEAD OF CAMERA

//right vector //x axis of camera local uvn tyoe system

const rightX = forward[1] * uz - forward[2] * uy;
const rightY = forward[2] * ux - forward[0] * uz;
const rightZ = forward[0] * uy - forward[1] * ux;
const rLength = Math.sqrt(rightX * rightX + rightY * rightY + rightZ * rightZ);
const right = [rightX / rLength, rightY / rLength, rightZ / rLength];  // normalised

// recompued up orthogonal uvn
const upX = right[1] * forward[2] - right[2] * forward[1];
const upY = right[2] * forward[0] - right[0] * forward[2];
const upZ = right[0] * forward[1] - right[1] * forward[0];
const upLength = Math.sqrt(upX * upX + upY * upY + upZ * upZ);
const upNormalized = [upX / upLength, upY / upLength, upZ / upLength];



/*
   const rotation = [
        [right[0], newUp[0], -forward[0], 0],
        [right[1], newUp[1], -forward[1], 0],
        [right[2], newUp[2], -forward[2], 0],
        [0, 0, 0, 1],
    ];
*/

/*

   const translation = [
        [1, 0, 0, -eye[0]],
        [0, 1, 0, -eye[1]],
        [0, 0, 1, -eye[2]],
        [0, 0, 0, 1],
    ];

*/


// view matrix  // rotation into translation
return[
    [right[0], upNormalized[0], -forward[0], 0],
    [right[1], upNormalized[1], -forward[1], 0],
    [right[2], upNormalized[2], -forward[2], 0],
    [
        -(right[0] * ex + right[1] * ey + right[2] * ez),
        -(upNormalized[0] * ex + upNormalized[1] * ey + upNormalized[2] * ez),
         forward[0] * ex + forward[1] * ey + forward[2] * ez,
            1,
    ]
  ];
}


function PerspectiveProjection(fovy, aspectRatio, znear, zfar) {
    const h = 1 / Math.tan((fovy * 0.5) * (Math.PI / 180)); // Convert fovy to radians
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
function createOrthographicMatrix(left, right, bottom, top, near, far) {
    const orthoMatrix = [
        2 / (right - left), 0, 0, -(right + left) / (right - left),
        0, 2 / (top - bottom), 0, -(top + bottom) / (top - bottom),
        0, 0, -2 / (far - near), -(far + near) / (far - near),
        0, 0, 0, 1
    ];

    return orthoMatrix;
}
*/

/*
transformAndProject(vertexbuffer, WorldMatrix, ViewMatrix, PerspectiveProjection){
    const transformedVertices = vertexbuffer.map(vertex => {
        const Homogeneous = [...vertex, 1];
        
        const combinedMatrix = MatrixMultiply(PerspectiveProjection, MatrixMultiply(ViewMatrix, WorldMatrix));
        let transformedVertex = MatrixVectorMultiply(combinedMatrix, Homogeneous);
        
        const w = transformedVertex[3]; if (Math.abs(w) > 0.0001) { 
            transformedVertex = transformedVertex.map(v => v / w); 
        }

        
        // Convert to screen space
        const screenX = Math.floor((transformedVertex[0] + 1) * width / 2);
        const screenY = Math.floor((1 - transformedVertex[1]) * height / 2);
        
        // Return in your exact format
        return [screenX, screenY, transformedVertex[2]]; 
    });

    // Return the complete array, exactly as in your code
    return transformedVertices;

    }*/

    function transformAndProject(vertexbuffer, WorldMatrix, ViewMatrix, PerspectiveProjection) {
        const transformedVertices = vertexbuffer.map(vertex => {
            const Homogeneous = [...vertex, 1];
            
            
            const combinedMatrix = MatrixMultiply(PerspectiveProjection, MatrixMultiply(ViewMatrix, WorldMatrix));
            let transformedVertex = MatrixVectorMultiply(combinedMatrix, Homogeneous);
            

            // this the problem 
            const w = transformedVertex[3]; if (Math.abs(w) > 0.0001) { 
                transformedVertex = transformedVertex.map(v => v / w); 
            }
    
        
            
            
            const screenX = Math.floor((transformedVertex[0] + 1) * width / 2);
            const screenY = Math.floor((1 - transformedVertex[1]) * height / 2);
            
          
            return [screenX, screenY, transformedVertex[2]];
        });
    
        return transformedVertices;
    }
    



/*
//the rule for nany line algo is , if m is 
less than 1 then always increment x and 
calculate y . if m more than 1 do opp incrment y and calculate x
*/



function drawlineH(x0,y0,x1,y1,r,g,b,a=255){
    if (x0>x1){   // the horizontal 4 octant
    [x0,x1] = [x1,x0];
    [y0,y1] = [y1,y0];
    }

    const dx = x1 - x0;
    const dy = y1 -y0;

     //covers slope in those 4 cocatnts //let dir = (dy < 0) ? -1 : 1;

    if(dy<0){
         dir = -1;
    }else{
        dir =1;
    }   
    
    if (dx!= 0){
        y = y0;
        let p = 2*dy -dx       // twice the half the pixel error for dx preload offset error // refer documentation
        for(let i=0; i<(dx +1); i++){
            setPixel(x0+i,y,r,g,b,a);                         //diagonal
            if(p>=0){
                y += dir;
                p = p -2*dx;
            }
            p = p+2*dy;
         }
    }
}


function drawlineV(x0,y0,x1,y1,r,g,b,a=255){
    if (y0>y1){
        [x0,x1] = [x1,x0];
        [y0,y1] = [y1,y0];
    }

    const dx = x1 - x0;
    const dy = y1 - y0;
    const dir = dx < 0 ? -1 : 1;
    
    if (dy !== 0){
        let x = x0;
        let p = 2*dx - dy;
        for(let i=0; i<(dy + 1); i++){
            setPixel(x, y0 + i, r,g,b,a);
            if(p>=0){
                x += dir;
                p = p - 2*dy;
            }
            p = p + 2*dx;
        }
    }
}




function DrawLine(x0,y0,x1,y1,r,g,b,a=255){
    if(Math.abs(x1 - x0) > Math.abs(y1-y0)){
        drawlineH(x0,y0,x1,y1,r,g,b,a)
    }else{
        drawlineV(x0,y0,x1,y1,r,g,b,a)
    }
}


/*

function drawTriangle(v0, v1, v2, r0, g0, b0, a0, r1, g1, b1, a1, r2, g2, b2, a2) {
    // Raster grid sort by y coords, scanline moves top to down
    if (v1[1] < v0[1]) {
        [v0, v1] = [v1, v0];
    }
    if (v2[1] < v1[1]) {
        [v1, v2] = [v2, v1];
    }
    if (v1[1] < v0[1]) {
        [v0, v1] = [v1, v0]; // v0 is the smallest
    }

    // Finding x at a particular y
    function interpolateY(x0, y0, x1, y1, y) {
        return x0 + ((x1 - x0) * (y - y0) / (y1 - y0));
    }

    function colorinterpolation(c0, c1, t) {
        return c0 + (c1 - c0) * t;
    }

    // Loop through the y values from the bottom vertex to the top vertex
    for (let y = Math.ceil(v0[1]); y <= Math.floor(v2[1]); y++) {
        let xleft, xright;
        let rleft, gleft, bleft, aleft;
        let rright, gright, bright, aright;

        if (y < v1[1]) {
            xleft = interpolateY(v0[0], v0[1], v1[0], v1[1], y);
            xright = interpolateY(v1[0], v1[1], v2[0], v2[1], y);

            rleft = colorinterpolation(r0, r1, (y - v0[1]) / (v1[1] - v0[1]));
            gleft = colorinterpolation(g0, g1, (y - v0[1]) / (v1[1] - v0[1]));
            bleft = colorinterpolation(b0, b1, (y - v0[1]) / (v1[1] - v0[1]));
            aleft = colorinterpolation(a0, a1, (y - v0[1]) / (v1[1] - v0[1]));

            rright = colorinterpolation(r0, r2, (y - v0[1]) / (v2[1] - v0[1]));
            gright = colorinterpolation(g0, g2, (y - v0[1]) / (v2[1] - v0[1]));
            bright = colorinterpolation(b0, b2, (y - v0[1]) / (v2[1] - v0[1]));
            aright = colorinterpolation(a0, a2, (y - v0[1]) / (v2[1] - v0[1]));
        } else {
            xleft = interpolateY(v1[0], v1[1], v2[0], v2[1], y);
            xright = interpolateY(v0[0], v0[1], v2[0], v2[1], y);

            rleft = colorinterpolation(r1, r2, (y - v1[1]) / (v2[1] - v1[1]));
            gleft = colorinterpolation(g1, g2, (y - v1[1]) / (v2[1] - v1[1]));
            bleft = colorinterpolation(b1, b2, (y - v1[1]) / (v2[1] - v1[1]));
            aleft = colorinterpolation(a1, a2, (y - v1[1]) / (v2[1] - v1[1]));

            rright = colorinterpolation(r0, r2, (y - v0[1]) / (v2[1] - v0[1]));
            gright = colorinterpolation(g0, g2, (y - v0[1]) / (v2[1] - v0[1]));
            bright = colorinterpolation(b0, b2, (y - v0[1]) / (v2[1] - v0[1]));
            aright = colorinterpolation(a0, a2, (y - v0[1]) / (v2[1] - v0[1]));
        }

        // Loop through x values between xleft and xright
        for (let x = Math.ceil(Math.min(xleft, xright)); x <= Math.floor(Math.max(xleft, xright)); x++) {
            let t = (x - Math.min(xleft, xright)) / (Math.max(xleft, xright) - Math.min(xleft, xright));

            let rvalue = colorinterpolation(rleft, rright, t);
            let gvalue = colorinterpolation(gleft, gright, t);
            let bvalue = colorinterpolation(bleft, bright, t);
            let avalue = colorinterpolation(aleft, aright, t);

            setPixel(x, y, rvalue, gvalue, bvalue, avalue); // Ensure setPixel is correctly implemented
        }
    }
}
*/


/*
function transformAndProject(transformedVertices, indexbuffer, colorbuffer) {
    for (let i = 0; i < indexbuffer.length; i++) {
        const [i0, i1, i2] = indexbuffer[i];
        
        // Direct access since vertices are already transformed
        // and in the correct order
        drawTriangle(
            transformedVertices[i0],
            transformedVertices[i1],
            transformedVertices[i2],
            ...colorbuffer[i0],  // spread the color components
            ...colorbuffer[i1],
            ...colorbuffer[i2]
        );
    }
}
*/



function interpolateY(x0, y0, x1, y1, y) {
    return x0 + (x1 - x0) * ((y - y0) / (y1 - y0));
}


function interpolateColor(c0, c1, t) {
    if (!c0 || !c1) {
        console.error("Invalid color values:", c0, c1);
        return [0, 0, 0]; 
    }
    return [
        Math.round(c0[0] + t * (c1[0] - c0[0])),
        Math.round(c0[1] + t * (c1[1] - c0[1])),
        Math.round(c0[2] + t * (c1[2] - c0[2])),
        Math.round(c0[3] + t * (c1[3] - c0[3])) /
    ];
}





function drawTriangle(v0, v1, v2, c0, c1, c2) {
    
    // top to botoom sort
    if (v1[1] < v0[1]) [v0, v1] = [v1, v0];
    if (v2[1] < v1[1]) [v1, v2] = [v2, v1];
    if (v1[1] < v0[1]) [v0, v1] = [v1, v0];

   
 

   
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

        /
        // color inteprolation
        for (let x = Math.ceil(Math.min(xLeft, xRight)); x <= Math.floor(Math.max(xLeft, xRight)); x++) {

            // Interpolating color along the horizontal scanline (between xLeft and xRight)
            const t = (x - Math.min(xLeft, xRight)) / (Math.max(xLeft, xRight) - Math.min(xLeft, xRight));
            const pixelColor = interpolateColor(cLeft, cRight, t);
            setPixel(x, y, pixelColor[0], pixelColor[1], pixelColor[2], pixelColor[3]);
        }
    }
}


/*

function drawanything(transformedVertices, indexbuffer, colorbuffer) {
    for (let i = 0; i < indexbuffer.length; i++) {
        const indices = indexbuffer[i];
        const [i0, i1, i2] = indices.slice(0, 3);  // Always use the first three indices for a triangle


        const c0 = colorbuffer[i0] || [0, 0, 0];  // Default to black if undefined
        const c1 = colorbuffer[i1] || [0, 0, 0];
        const c2 = colorbuffer[i2] || [0, 0, 0];

        console.log(`Using colors: c0 = ${c0}, c1 = ${c1}, c2 = ${c2}`);
        
        if (!c0 || !c1 || !c2) {
            console.error("Invalid color values detected:", c0, c1, c2);
        }

        // Draw the first triangle (always 3 vertices in the first pass)
        drawTriangle(
            transformedVertices[i0],
            transformedVertices[i1],
            transformedVertices[i2],
            colorbuffer[i0],  // Spread the color components for vertex i0
            colorbuffer[i1],  // Spread the color components for vertex i1
            colorbuffer[i2]   // Spread the color components for vertex i2
        );


        // If it's a quad (i.e., 4 indices), draw the second triangle
        if (indices.length === 4) {
            const [i3] = indices.slice(3, 4);  // Get the 4th index for the second triangle
            drawTriangle(
                transformedVertices[i0],
                transformedVertices[i2],
                transformedVertices[i3],
                colorbuffer[i0],
                colorbuffer[i2],
                colorbuffer[i3]
            );

        }

        console.log("Colorbuffer values:", colorbuffer[i0], colorbuffer[i1], colorbuffer[i2]);
    }
}
*/






function drawanything(transformedVertices, indexbuffer, colorbuffer) {
    for (let i = 0; i < indexbuffer.length; i++) {
        const indices = indexbuffer[i];
        const [i0, i1, i2] = indices.slice(0, 3);
        

        console.log("Indexbuffer values:", indices);

        // agar mer epe undefined aayege papieline  me default to balck
        const c0 = colorbuffer[i0] || [0, 0, 0];
        const c1 = colorbuffer[i1] || [0, 0, 0];
        const c2 = colorbuffer[i2] || [0, 0, 0];
        const c3 = colorbuffer[i3] || [0, 0, 0];
        
        console.log(`Using colors: c0 = ${c0}, c1 = ${c1}, c2 = ${c2}, c3 = ${c3}`);          
        if (!c0 || !c1 || !c2 || !c3) {
            console.error("Invalid color values detected:", c0, c1, c2);
        }

        drawTriangle(
            transformedVertices[i0],
            transformedVertices[i1],
            transformedVertices[i2],
            c0, 
            c1,  
            c2   
        );

       
        // triangulation kinda
        if (indices.length === 4) {
            const [i3] = indices.slice(3, 4);
            
            console.log(" colors for second triangle:", colorbuffer[i0], colorbuffer[i2], colorbuffer[i3]);
           
            drawTriangle(
                transformedVertices[i0],
                transformedVertices[i2],
                transformedVertices[i3],
                c0,  
                c2, 
                c3,
            );
        }

        console.log("Colorbuffer values:", colorbuffer[i0], colorbuffer[i1], colorbuffer[i2]);
    }
}



const vertexbuffer = [
    [-1, -1, -1],  
    [1, -1, -1],   
    [1, 1, -1],    
    [-1, 1, -1],   
    [-1, -1, 1],   
    [1, -1, 1],    
    [1, 1, 1],     
    [-1, 1, 1], 
    
    [2, -2, -2],   
    [4, -2, -2],   
    [3, 1, -2],
];

const colorbuffer = [
    [255, 0, 0],    
    [0, 255, 0],   
    [0, 0, 255],    
    [255, 255, 0],  
    [255, 0, 255],  
    [0, 255, 255],  
    [255, 165, 0],  
    [128, 0, 128],
    
    [255, 255, 255], 
    [255, 255, 255], 
    [255, 255, 255],
    
];


const indexbuffer = [
    [0, 1, 2], [0, 2, 3], 
    [4, 5, 6], [4, 6, 7], 
    [0, 1, 5], [0, 5, 4], 
    [1, 2, 6], [1, 6, 5], 
    [2, 3, 7], [2, 7, 6], 
    [3, 0, 4], [3, 4, 7],
    
    [8,9,10],
];





/*
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
*/

/*
function render() {
    
    ClearFrameBuffer(0, 0, 0);

    const eye = [0, 0, 3];
    const center = [0, 0, 0];
    const up = [0, 1, 0];

    const worldMat = WorldMatrix(0, 0, 0, 0, angle, 0, 1, 1, 1);
    const viewMat = ViewMatrix(eye, center, up);
    const projMat = PerspectiveProjection(45, canvas.width / canvas.height, 0.1, 1000);

    const transformedVertices = transformedrender(vertexbuffer, worldMat, viewMat, projMat);

      
       
       // Draw vertices as larger points first
       transformedVertices.forEach((vertex, i) => {
           const color = colorBuffer[i];
           // Draw larger points (5x5 pixels) at each vertex
           for(let y = -5; y <= 5; y++) {
               for(let x = -5; x <= 5; x++) {
                   setPixel(vertex[0] + x, vertex[1] + y, color[0], color[1], color[2]);
               }
           }
       });
       
       // Draw lines between vertices to see triangle outline
       const indices = indexBuffer[0];
       for(let i = 0; i < indices.length; i++) {
           const start = transformedVertices[indices[i]];
           const end = transformedVertices[indices[(i + 1) % 3]];
           DrawLine(start[0], start[1], end[0], end[1], 255, 255, 255);
       }


    Trianglesrender(transformedVertices, indexBuffer, colorBuffer);

    ctx.putImageData(framebuffer, 0, 0);

    angle += 1;
    requestAnimationFrame(render);
    
}
*/

/*

function render() {
    clearFrameBuffer(0, 0, 0, 255);

    const eye = [0, 0, 3];
    const center = [0, 0, 0];
    const up = [0, 1, 0];
    const aspectRatio = width / height;
    const projection = PerspectiveProjection(90, aspectRatio, 0.1, 100);

    const translation = translationMatrix(0, 0, -2);
    const rotation = RotationMatrixY(angle);

    const model = MatrixMultiply( MatrixMultiply(translation, rotation));
    const worldmat = MatrixMultiply(model,ViewMatrix);
    const finalmat = MatrixMultiply(worldmat ,PerspectiveProjection)

    const transformedVertices = vertexbuffer.map(vertex =>
        MatrixVectorMultiply(finalmat, [...vertex, 1])
    ).map(v => [
        Math.floor((v[0] / v[3] + 1) * width / 2),
        Math.floor((1 - v[1] / v[3]) * height / 2),
    ]);

    // Pass the individual colors for each vertex
    drawTriangle(
        transformedVertices[0], transformedVertices[1], transformedVertices[2],colorBuffer[0],colorBuffer[1],colorBuffer[2]  
    );
    


    ctx.putImageData(framebuffer, 0, 0);

    angle += 1;
    requestAnimationFrame(render);
}

let angle = 0;
render();
*/



//render , draw trainglkes adn drwa me 


/*
function render(){

    ClearCanvas();
    clearFrameBuffer(0, 0, 0);  
    

    
    const worldMatrix = worldMatrix(0, 0, 0, 0, 0, 0, 1, 1, 1);  
    const viewMatrix = viewMatrix([0, 0, -5], [0, 0, 0], [0, 1, 0]); 
    const projectionMatrix = PerspectiveProjection(90, width / height, 0.1, 100);

    // Transform and project the vertices
    const transformedVertices = transformAndProject(vertexbuffer, worldMatrix, viewMatrix, projectionMatrix);
    
    for (let i = 0; i < indexbuffer.length; i++) {
        const indices = indexbuffer[i];
        const v0 = transformedVertices[indices[0]];
        const v1 = transformedVertices[indices[1]];
        const v2 = transformedVertices[indices[2]];

        // Here you can draw the triangle with interpolated colors
        drawTriangle(v0, v1, v2, ...colorbuffer[0], ...colorbuffer[1], ...colorbuffer[2]);
    }

    
    ctx.putImageData(framebuffer, 0, 0);
    requestAnimationFrame(render)
}   

render();

*/

let angle = 0;

function render() {
    clearFrameBuffer(0, 0, 0, 255);

    const aspectRatio = width / height;
    const projection = PerspectiveProjection(90, aspectRatio, 0.1, 1000);

    const translation = translationMatrix(0, 0, 5);
    const rotation = RotationMatrixY(angle);
    const transform = MatrixMultiply(MatrixMultiply(projection, translation), rotation);

    const transformedVertices = vertexbuffer.map(vertex => {

        console.log('Before transformation:', vertex);
        const v = MatrixVectorMultiply(transform, [...vertex, 1]);
        console.log('After transformation:', v);
        const w = v[3];
        return [
            Math.floor((v[0] / w) * width / 2 + width / 2),
            Math.floor((-v[1] / w) * height / 2 + height / 2),
            v[2] / w,  // depth to handle
        ];
    });
    console.log("Transformed Vertices:", transformedVertices);

    drawanything(transformedVertices, indexbuffer, colorbuffer);

    ctx.putImageData(framebuffer, 0, 0);

    angle += 1;
    requestAnimationFrame(render);
}

render();

  // shaktiman chalkra


/*
let speed = 0.1; 

function render() {
    
    ClearCanvas();
    clearFrameBuffer(0, 0, 0);

    
    console.log("Vertex Buffer Before Render:", vertexbuffer);
    console.log("Index Buffer Before Render:", indexbuffer);

    //const worldMatrix = WorldMatrix(0, 0, 0, 0, 0, 0, 1, 1, 1);
    rotationangle += speed;
    speed += 0.1;
    const rotationmatrix = RotationMatrixZ(rotationangle);
    const viewMatrix = ViewMatrix([0, 0, -5], [0, 0, 0], [0, 1, 0]);
    const projectionMatrix = PerspectiveProjection(90, width / height, 0.1, 100);

    console.log("World Matrix:", rotationmatrix);
    console.log("View Matrix:", viewMatrix);
    console.log("Projection Matrix:", projectionMatrix);
    
    const transformedVertices = transformAndProject(vertexbuffer, rotationmatrix, viewMatrix, projectionMatrix);


    console.log("Transformed Vertices After Projection:", transformedVertices)

    // Draw all triangles/quads in the scene
    drawanything(transformedVertices, indexbuffer, colorbuffer);

    
    ctx.putImageData(framebuffer, 0, 0);

    rotationangle += 1;
    requestAnimationFrame(render);
}


let rotationangle = 0;
render();
*/



// have to check camera far positions

/*
//let speed =0.1;
function render() {
    
    ClearCanvas();
    clearFrameBuffer(0, 0, 0);

    

    //const worldMatrix = WorldMatrix(0, 0, 0, 0, 0, 0, 1, 1, 1);
    //rotationangle += speed;
    //speed += 0.1;
    const rotationmatrix = RotationMatrixZ(rotationangle);
    const viewMatrix = ViewMatrix([0, 0, -5], [0, 0, 0], [0, 1, 0]);
    const projectionMatrix = PerspectiveProjection(90, width / height, 1, 100);

    
    const transformedVertices = transformAndProject(vertexbuffer, rotationmatrix, viewMatrix, projectionMatrix);

    // Draw all triangles/quads in the scene
    drawanything(transformedVertices, indexbuffer, colorbuffer);

    
    ctx.putImageData(framebuffer, 0, 0);

    rotationangle +=1;
    
    requestAnimationFrame(render);
}


let rotationangle = 0;
render();
*/

/*

function render() {
    // Clear the canvas at the start of each frame
    ClearCanvas();
    clearFrameBuffer(0, 0, 0);  // Optional: Clear framebuffer

    // Set the transformation matrices for the world, view, and projection
    const worldMatrix = WorldMatrix(0, 0, 0, 0, 0, 0, 1, 1, 1);  // Example: Identity transformation
    const viewMatrix = ViewMatrix([0, 0, -5], [0, 0, 0], [0, 1, 0]);  // Camera setup
    const projectionMatrix = PerspectiveProjection(90, width / height, 0.1, 100);

    // Transform and project the vertices
    const transformedVertices = transformAndProject(vertexbuffer, worldMatrix, viewMatrix, projectionMatrix);

    // Loop through the indexBuffer to draw each triangle
    for (let i = 0; i < indexbuffer.length; i++) {
        const indices = indexbuffer[i];
        const v0 = transformedVertices[indices[0]];
        const v1 = transformedVertices[indices[1]];
        const v2 = transformedVertices[indices[2]];

        // Draw the first triangle
        drawTriangle(v0, v1, v2, colorbuffer[0], colorbuffer[1], colorbuffer[2]);

        // If there's a fourth vertex for the second triangle
        if (indices.length === 4) {
            const v3 = transformedVertices[indices[3]];
            // Draw the second triangle
            drawTriangle(v0, v2, v3, colorbuffer[0], colorbuffer[2], colorbuffer[3]);
        }
    }

    // Update the canvas framebuffer
    ctx.putImageData(framebuffer, 0, 0);
}

render();
*/


/*

function render() {
    clearFrameBuffer(0, 0, 0, 255);

    const aspectRatio = width / height;
    const projection = PerspectiveProjection(270, aspectRatio, 0.1, 100);

    const translation = translationMatrix(0, 0, -2);
    const rotation = RotationMatrixY(angle);
    const transform = MatrixMultiply(projection, MatrixMultiply(translation, rotation));

    const transformedVertices = vertexbuffer.map(vertex =>
        MatrixVectorMultiply(transform, [...vertex, 1])
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



/*
let rotationangle = 0;
function render() {
    
    ClearCanvas();
    clearFrameBuffer(0, 0, 0);

    

    //const worldMatrix = WorldMatrix(0, 0, 0, 0, 0, 0, 1, 1, 1);
    rotationangle += speed;
    speed += 0.1;
    const rotationmatrix = RotationMatrixZ(rotationangle);
    const eye = [0,0,-1000];
    const center = [100000 , -5333 ,0];
    const viewMatrix = ViewMatrix(eye, center, [0, 1, 0]);
    console.log("Camera Position:", eye);
    console.log("Camera Position center:", center);

    const projectionMatrix = PerspectiveProjection(60, width / height, 0.1, 100);

    
    const transformedVertices = transformAndProject(vertexbuffer, worldMatrix, viewMatrix, projectionMatrix);

    // Draw all triangles/quads in the scene
    drawanything(transformedVertices, indexbuffer, colorbuffer);

    
    ctx.putImageData(framebuffer, 0, 0);

    //rotationangle += 1;

    requestAnimationFrame(render);
}



render();
*/


/*
function render() {
    ClearCanvas();
    clearFrameBuffer(0, 0, 0);  // Black background

    const eye = [0, 0, 3];      // Camera position
    const center = [0, 0, 0];   // Looking at center
    const up = [0, 1, 0];       // Up direction

    // Create transformation matrices
    const worldMat = WorldMatrix(0, 0, 0, angle, angle * 0.5, 0, 1, 1, 1);
    const viewMat = ViewMatrix(eye, center, up);
    const projMat = PerspectiveProjection(45, canvas.width / canvas.height, 0.1, 1000);

    // Transform vertices
    const transformedVertices = transformAndProject(vertexbuffer, worldMat, viewMat, projMat);

    // Draw vertices as points
    transformedVertices.forEach(vertex => {
        // Draw a small point for each vertex (3x3 pixels)
        for(let y = -1; y <= 1; y++) {
            for(let x = -1; x <= 1; x++) {
                setPixel(
                    Math.floor(vertex[0] + x),
                    Math.floor(vertex[1] + y),
                    255, 0, 0  // Red vertices
                );
            }
        }
    });

    // Draw edges
    edgebuffer.forEach(([startIdx, endIdx]) => {
        const start = transformedVertices[startIdx];
        const end = transformedVertices[endIdx];
        
        DrawLine(
            Math.floor(start[0]), Math.floor(start[1]),
            Math.floor(end[0]), Math.floor(end[1]),
            255, 255, 255  // White edges
        );
    });

    // Update canvas
    ctx.putImageData(framebuffer, 0, 0);

    angle += 0.01;  // Rotate cube
    requestAnimationFrame(render);
}

// Start the rendering
let angle = 0;
render();

*/