const key={
    w:false,
    a:false,
    s:false,
    d:false
}
window.addEventListener("keydown",e=>{
    switch(e.code){
        case "KeyW":
            key.w=true;
            break;
        case "KeyA":
            key.a=true;
            break;
        case "KeyS":
            key.s=true;
            break;
        case "KeyD":
            key.d=true;
            break;
        case "KeyP":
            projectionId=(projectionId+1)%projectionName.length;
            break;
    }
});
window.addEventListener("keyup",e=>{
    switch(e.code){
        case "KeyW":
            key.w=false;
            break;
        case "KeyA":
            key.a=false;
            break;
        case "KeyS":
            key.s=false;
            break;
        case "KeyD":
            key.d=false;
            break;
    }
});
function keycontrol(){
    let v=[0,0];
    if(key.w){
        v[1]+=1;
    }
    if(key.a){
        v[0]+=1;
    }
    if(key.s){
        v[1]-=1;
    }
    if(key.d){
        v[0]-=1;
    }
    const s=20*Math.hypot(v[0],v[1]);
    if(s>0){
        v[0]*=1/s;
        v[1]*=-1/s;
    }
    moveVector=v;
}
render();