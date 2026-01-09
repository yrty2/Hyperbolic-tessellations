let tiling=new Float32Array([7,3]);
let moveVector=[0,0];
const canvas=document.querySelector(".canvas");
canvas.width=screen.width;
canvas.height=screen.height;
const ctx=canvas.getContext("2d");
const points=[];
let polygon=[];
const tiles=[];
let line=[];
const projectionName=["ポアンカレの円板モデル","クラインの円板モデル","ガンズモデル","上半平面モデル"];
let projectionId=0;
const center=new hyperbolic(0,0);
ctx.font="30px sans-serif";
function render(){
    ctx.fillStyle="#000000";
    if(projectionId==3 || projectionId==2){
        ctx.fillRect(0,0,canvas.width,canvas.height);
    }else{
    ctx.clearRect(0,0,canvas.width,canvas.height);
    ctx.beginPath();
    ctx.arc(canvas.width/2,canvas.height/2,canvas.height/2,0,2*Math.PI);
    ctx.fill();
    ctx.closePath();
    }
    for(const p of points){
        p.translate(moveVector[0],moveVector[1]);
        ctx.fillStyle=`hsl(${math.mod(180*p.arg/Math.PI,360)},100%,20%)`;
        plot(p);
    }
    for(const l of line){
        l.center.translate(moveVector[0],moveVector[1]);
        for(const q of l.pos){
            q.translate(moveVector[0],moveVector[1]);
        }
    }
    for(const l of line){
        if(projectionId!=3 || l.center.y<0.7){
        lineList(l);
        }
    }
    ctx.fillStyle="#000000";
    ctx.fillText(`{${tiling[0]},${tiling[1]}}`,30,130);
    ctx.fillText(`投影:${projectionName[projectionId]}`,30,90);
    ctx.fillText(`WASDキーで移動`,30,canvas.height-90);
    ctx.fillText(`Pキーで投影を変更`,30,canvas.height-50);
    keycontrol();
    requestAnimationFrame(render);
}
function tessellation(n,m,iter){
    let center=[];
    let prec=[c32.zero()];
    let s=2*Math.atanh(Math.sqrt((1/Math.tan(Math.PI/m)-Math.tan(Math.PI/n))/(1/Math.tan(Math.PI/m)+Math.tan(Math.PI/n))));
    tiles.push(c32.zero());
    ngon(n,s,0,0);
    const p=c32.zero();
    const H=polygon[polygon.length-1].pos.slice();
    let polygonList=[H];
    function create(level){
        const stockpile=[];
        const stockp=[];
    for(let i=0; i<polygonList.length; ++i){
    for(let k=0; k<n; ++k){
        const prev=polygonList[i];
        const cent=hyperbolicGeometry.reflection(prev[k].z,prev[(k+1)%n].z,prec[i]);
        let hantei=true;
        for(const c of [...center,...stockp,...prec]){
            if(hyperbolicGeometry.distance(cent,c)<0.01){
                hantei=false;
                break;
            }
        }
        if(hantei){
            stockp.push(cent);
            const poly=hyperbolicGeometry.reflectionPolygon(prev[k].z,prev[(k+1)%n].z,prev);
            stockpile.push(poly);
            polygon.push({center:cent,color:`hsl(${60*level},45%,60%)`,pos:poly});
        }
    }
        }
        center.push(...prec);
        prec=stockp;
        polygonList=stockpile;
    }
    for(let k=1; k<=iter; ++k){
        create(k);
    }
    instantiate();
}
tessellation(tiling[0],tiling[1],5);
instantiate();
function instantiate(){
    line=[];
    for(const p of polygon){
        line.push({center:new hyperbolic(p.center),color:p.color,pos:[]});
        for(let k=0; k<p.pos.length; ++k){
            const geo=p.pos[k].geodesic(p.pos[(k+1)%p.pos.length]);
            for(const g of geo){
            line[line.length-1].pos.push(new hyperbolic(g));
            }
        }
    }
}
//generate(50);
function generate(d){
    for(let x=-d; x<=d; ++x){
    for(let y=-d; y<=d; ++y){
    points.push(new hyperbolic(x/100,y/100));
    }}
}
function ngon(n,s,theta,color,...z){
    polygon.push({center:c32.zero(),color:`hsl(${color},45%,60%)`,pos:[]});
    for(let k=0; k<n; ++k){
        const h=new hyperbolic(s,0);
        h.R(2*k/n*Math.PI);
        if(theta){
        h.R(theta);
        }
        for(let k=0; k<z.length; ++k){
            h.translate(z[k][0],z[k][1]);
        }
        polygon[polygon.length-1].pos.push(h);
    }
}
function plot(p){
    p.plot(canvas,ctx);
}
function lineList(linedata){
        ctx.beginPath();
        for(const p of linedata.pos){
        let q=p.z;
        switch (projectionId){
            case 1:
                q=p.klein;
                break;
            case 2:
                q=c32.prod(p.disk(128),32);
                break;
            case 3:
                q=p.upperhalf;
                break;
        }
        ctx.lineTo((canvas.height*q[0]+canvas.width)/2,canvas.height*(-q[1]+1)/2);
        }
        ctx.fillStyle=linedata.color;
        ctx.fill();
        ctx.stroke();
        ctx.closePath();
}
function change(p,q,iter){
    if((p-2)*(q-2)>4){
    tiling=[p,q];
    polygon=[];
    tessellation(tiling[0],tiling[1],iter);
    }else{
        if((p-2)*(q-2)==4){
            alert("それはユークリッド幾何学です！");
        }else{
            if(p>1 && q>1){
                alert("それは球面幾何学です！");
            }else{
                alert("値が不適切です！");
            }
        }
    }
}