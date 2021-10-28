var kage_export = (function (exports) {
  'use strict';

  class Buhin{
    constructor(number){
      this.hash = {};
    }
    set(name, data){ // void
      this.hash[name] = data;
    }
    push(name, data){ // void
      this.set(name, data);
    }search(name){ // string
      if(this.hash[name]){
        return this.hash[name];
      }
      return ""; // no data
    }
  }

  class Polygons{
    // method
    constructor(){
      this.array = new Array();
    }
   	clear(){ // void
      this.array = new Array();
    }
    push(polygon){ // void
      // only a simple check
      var minx = 200;
      var maxx = 0;
      var miny = 200;
      var maxy = 0;
      var error = 0;
      for(var i = 0; i < polygon.array.length; i++){
        if(polygon.array[i].x < minx){
          minx = polygon.array[i].x;
        }
        if(polygon.array[i].x > maxx){
          maxx = polygon.array[i].x;
        }
        if(polygon.array[i].y < miny){
          miny = polygon.array[i].y;
        }
        if(polygon.array[i].y > maxy){
          maxy = polygon.array[i].y;
        }
        if(isNaN(polygon.array[i].x) || isNaN(polygon.array[i].y)){
          error++;
        }
      }
      if(error == 0 && minx != maxx && miny != maxy && polygon.array.length >= 3){
        var newArray = new Array();
        newArray.push(polygon.array.shift());
        while(polygon.array.length != 0){
          var temp = polygon.array.shift();
          //if(newArray[newArray.length - 1].x != temp.x ||
          //   newArray[newArray.length - 1].y != temp.y){
            newArray.push(temp);
          //}
        }
        if(newArray.length >= 3){
          polygon.array = newArray;
          this.array.push(polygon);
        }
      }
    }
    concat(polygons){
      for(let polygon of polygons.array){
        this.push(polygon);
      }
    }
    generateSVG(){ // string
      var buffer = "";
      buffer += "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" version=\"1.1\" baseProfile=\"full\" viewBox=\"0 0 200 200\" width=\"200\" height=\"200\">\n";
      for(var i = 0; i < this.array.length; i++){
        buffer += "<path d=\"";
        buffer += this.array[i].get_sub_path_svg();
        buffer += "\" fill=\"black\" />\n";
      }
      buffer += "</svg>\n";
      return buffer;
    }
    generateSVG2(){ // string
      var buffer = "";
      buffer += "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" version=\"1.1\" baseProfile=\"full\" viewBox=\"0 0 200 200\" width=\"200\" height=\"200\">\n";
      buffer += "<path d=\"";
      buffer += this.get_path_svg();
      buffer += "\" fill=\"black\" />\n";
      buffer += "</svg>\n";
      return buffer;
    }
    get_path_svg(){
      var buffer = "";
      for(var i = 0; i < this.array.length; i++){
        buffer += this.array[i].get_sub_path_svg();
      }
      return buffer;
    }
    get_path_svg_font(){
      var buffer = "";
      for(var i = 0; i < this.array.length; i++){
        buffer += this.array[i].get_sub_path_svg_font();
      }
      return buffer;
    }
    generateEPS(){ // string
      var buffer = "";
      buffer += "%!PS-Adobe-3.0 EPSF-3.0\n";
      buffer += "%%BoundingBox: 0 -208 1024 816\n";
      buffer += "%%Pages: 0\n";
      buffer += "%%Title: Kanji glyph\n";
      buffer += "%%Creator: GlyphWiki powered by KAGE system\n";
      buffer += "%%CreationDate: " + new Date() + "\n";
      buffer += "%%EndComments\n";
      buffer += "%%EndProlog\n";
      
      for(var i = 0; i < this.array.length; i++){
        for(var j = 0; j < this.array[i].array.length; j++){
          buffer += (this.array[i].array[j].x * 5) + " " + (1000 - this.array[i].array[j].y * 5 - 200) + " ";
          if(j == 0){
            buffer += "newpath\nmoveto\n";
          } else {
            buffer += "lineto\n";
          }
        }
        buffer += "closepath\nfill\n";
      }
      buffer += "%%EOF\n";
      return buffer;
    }

  }

  class Polygon{
    // resolution : 0.0001
    constructor(number){
      this.array = new Array();
      // initialize
      if(number){
        for(var i = 0; i < number; i++){
          this.push(0, 0, 0);
        }
      }
    }
    push(x, y, off){ // void
      var temp = new Object();
      temp.x = Math.round(x*10000)/10000;
      temp.y = Math.round(y*10000)/10000;
      if(off != 1 && off != 2){
        off = 0;
      }
      temp.off = off;
      this.array.push(temp);
    }
    push2(p, off){
      let [x, y] = p;
      this.push(x, y, off);
    }
    set(index, x, y, off){ // void
      this.array[index].x = Math.round(x*10000)/10000;
      this.array[index].y = Math.round(y*10000)/10000;
      if(off != 1 && off != 2){
        off = 0;
      }
      this.array[index].off = off;
    }
    to_font1000(){
      for(var j = 0; j < this.array.length; j++){
        this.array[j].x = this.array[j].x*5;
        this.array[j].y = this.array[j].y*-5 - 200;
      }
    }
    reverse(){ // void
      this.array.reverse();
    }
    concat(poly){ // void
      this.array = this.array.concat(poly.array);
    }
    shift(){ // void
      this.array.shift();
    }
    unshift(x, y, off){ // void
      var temp = new Object();
      temp.x = Math.round(x*10000)/10000;
      temp.y = Math.round(y*10000)/10000;
      if(off != 1 && off != 2){
        off = 0;
      }
      temp.off = off;
      this.array.unshift(temp);
    }
    get_sub_path_svg(){
      let buffer = "";
      buffer += "M";
      buffer += this.array[0].x + "," + this.array[0].y + " ";
      let mode = "";
      for(var j = 1; j < this.array.length; j++){
        if(this.array[j].off == 1 && mode != "Q"){
          buffer += "Q"; mode = "Q";
        } else if(this.array[j-1].off == 0 && this.array[j].off == 2 && mode != "C"){
          buffer += "C"; mode = "C";
        } else if(this.array[j-1].off == 0 && this.array[j].off == 0 && mode != "L"){
          buffer += "L"; mode = "L";
        }
        buffer += this.array[j].x + "," + this.array[j].y + " ";
      }
      buffer += "Z";
      return buffer;
    }
    get_sub_path_svg_font(){
      let buffer = "";
      buffer += "M";
      buffer += (this.array[0].x*5) + "," + (800-this.array[0].y*5) + " ";
      let mode = "";
      for(var j = 1; j < this.array.length; j++){
        if(this.array[j].off == 1 && mode != "Q"){
          buffer += "Q"; mode = "Q";
        } else if(this.array[j-1].off == 0 && this.array[j].off == 2 && mode != "C"){
          buffer += "C"; mode = "C";
        } else if(this.array[j-1].off == 0 && this.array[j].off == 0 && mode != "L"){
          buffer += "L"; mode = "L";
        }
        buffer += (this.array[j].x*5) + "," + (800-this.array[j].y*5) + " ";
      }
      buffer += "Z";
      return buffer;
    }
  }

  const bez_cir = 4*(Math.sqrt(2)-1)/3;
  //a constant for drawing circles with Bezier curves

  //width functions (using circle)
  function widfun(t, x1, y1, x2, y2, wid){
    const len = Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    const p = 1 + (100/len);
    return (  (Math.sqrt(p*p+(p-1)*(p-1)-(p-t)*(p-t))-(p-1))*0.778+0.222  )*wid;
  }

  function widfun_d(t, x1, y1, x2, y2, wid){
    const len = Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    const p = 1 + (100/len);
    return wid*0.778*0.5*2*(p-t) / Math.sqrt(p*p+(p-1)*(p-1)-(p-t)*(p-t));
  }

  function widfun_stop(t, x1, y1, x2, y2, wid){
    const len = Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    const p = 1 + (100/len);
    return (  (Math.sqrt(p*p+(p-1)*(p-1)-(p-t)*(p-t))-(p-1))*0.878+0.222  )*wid;
  }

  function widfun_stop_d(t, x1, y1, x2, y2, wid){
    const len = Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    const p = 1 + (100/len);
    return wid*0.878*0.5*2*(p-t) / Math.sqrt(p*p+(p-1)*(p-1)-(p-t)*(p-t));
  }

  //fat version (used in cubic bezier)
  function widfun_fat(t, x1, y1, x2, y2, wid){
    const len = Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    const p = 1+ (40/len);
    return (  (Math.sqrt(p*p + (p-1)*(p-1) - (p-t)*(p-t)) - (p-1)  )*0.778+0.222  )*wid;
  }

  function widfun_fat_d(t, x1, y1, x2, y2, wid){
    const len = Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    const p = 1+ (40/len);
    return wid*0.778*0.5*2*(p-t) / Math.sqrt(p*p + (p-1)*(p-1) - (p-t)*(p-t));
  }

  function get_dir(x, y){
    if (y==0){
      if(x<0){
        return {sin: 0, cos: -1};
      }else {
        return {sin: 0, cos:  1};
      }
    }else if(x==0){
      if(y<0){
        return {sin: -1, cos: 0};
      }else {
        return {sin:  1, cos: 0};
      }
    }else {
      const rad = Math.atan2(y, x);
      return {sin:  Math.sin(rad), cos: Math.cos(rad)};
    }
  }
  const DIR_POSX = {sin: 0, cos: 1};
  const DIR_NEGX = {sin: 0, cos: -1};

  function moved_point(x, y, dir, delta){
    return [x + delta*dir.cos, y + delta*dir.sin];
  }
  function get_extended_dest(destX, destY, srcX, srcY, delta) {
    const dir = get_dir(destX - srcX, destY - srcY);
    return moved_point(destX, destY, dir, delta);
  }

  function get_extended_dest_wrong(destX, destY, srcX, srcY, delta) {
    //The process for lines directed exactly in the negative x-direction or y-direction is not correct, so it's named as "wrong".
    var destX_new = destX;
    var destY_new = destY;
    if (srcX == destX) {
      destY_new = destY + delta;
    }
    else if (srcY == destY) {
      destX_new = destX + delta;
    }
    else {
      var v;
      const rad = Math.atan((destY - srcY) / (destX - srcX));
      if (srcX < destX) { v = 1; } else { v = -1; }
      destX_new = destX + delta * Math.cos(rad) * v;
      destY_new = destY + delta * Math.sin(rad) * v;
    }
    return [destX_new, destY_new]
  }

  function unit_normal_vector(ix, iy) {//to the right(clockwise (in the display coordinate))
    var ia, ib;
    // line SUICHOKU by vector
    if (ix != 0 && iy != 0) {
      const ir = Math.atan(iy / ix * -1.0);
      ia = Math.sin(ir);
      ib = Math.cos(ir);
    }
    else if (ix == 0) {
      if (iy < 0) {
        ia = -1;
      } else {
        ia = 1;
      }
      ib = 0;
    }
    else {
      ia = 0;
      ib = 1;
    }
    //reverse if vector is going 2nd/3rd quadrants
    if (ix <= 0) {
      ia = ia * -1;
      ib = ib * -1;
    }
    return [ia, ib];
  }

  function vector_to_len(v, l){
    const len=Math.sqrt(v[0]*v[0]+v[1]*v[1]);
    return [v[0]*l/len,v[1]*l/len];
  }

  function get_rad(x, y) {
    var rad;
    if (x == 0) {
      if (y > 0) {
        rad = Math.PI / 2;
      } else {
        rad = -Math.PI / 2;
      }
    } else {
      rad = Math.atan(y / x);
      if (x < 0) { rad += Math.PI; }
    }
    return rad;
  }

  function rad_to_vector(rad) {
    return [Math.cos(rad), Math.sin(rad)];
  }

  function stretch_bezier_end(bez, t){
    const start = [bez[0][0],
                   bez[0][1]];
    const c1 = [(1-t) * bez[0][0]+t * bez[1][0],
                (1-t) * bez[0][1]+t * bez[1][1]];
    const c2 = [(1-t) * (1-t) * bez[0][0] + 2.0 * t * (1-t) * bez[1][0] + t * t * bez[2][0],
                (1-t) * (1-t) * bez[0][1] + 2.0 * t * (1-t) * bez[1][1] + t * t * bez[2][1]];
    const end = [(1-t) * (1-t) * (1-t) * bez[0][0] + 3 * t * (1-t) * (1-t) * bez[1][0] + 3 * t * t * (1-t) * bez[2][0] + t * t * t * bez[3][0],
                 (1-t) * (1-t) * (1-t) * bez[0][1] + 3 * t * (1-t) * (1-t) * bez[1][1] + 3 * t * t * (1-t) * bez[2][1] + t * t * t * bez[3][1],];
    return [start, c1, c2, end];
  }
  function bezier_to_y(bez, y){
    const res = shorten_bezier_to_y(bez, y);
    if(res){return res;}else {
      return extend_bezier_to_y(bez, y);
    }
  }
  function extend_bezier_to_y(bez, y) {
    const a =     bez[3][1] - 3 * bez[2][1] + 3 * bez[1][1] - bez[0][1];
    const b = 3 * bez[2][1] - 6 * bez[1][1] + 3 * bez[0][1];
    const c = 3 * bez[1][1] - 3 * bez[0][1];
    const d =     bez[0][1];
    const yy = solveCubic(a, b, c, d - y);
    yy.sort(function (a, b) {//ascending order
      return a - b;
    });
    for (let i of yy) {
      if (i > 1) {
        return stretch_bezier_end(bez, i);
      }
    }
    return false;
  }

  function shorten_bezier_to_y(bez, y) {
    const a =     bez[3][1] - 3 * bez[2][1] + 3 * bez[1][1] - bez[0][1];
    const b = 3 * bez[2][1] - 6 * bez[1][1] + 3 * bez[0][1];
    const c = 3 * bez[1][1] - 3 * bez[0][1];
    const d =     bez[0][1];
    const yy = solveCubic(a, b, c, d - y);
    yy.sort(function (a, b) {//descending order
      return b - a;
    });
    for (let i of yy) {
      if (0 < i && i < 1) {
        return stretch_bezier_end(bez, i);
      }
    }
    return false;
  }

  function solveCubic(a, b, c, d) {
    if (Math.abs(a) < 1e-8) { // Quadratic case, ax^2+bx+c=0
        a = b; b = c; c = d;
        if (Math.abs(a) < 1e-8) { // Linear case, ax+b=0
            a = b; b = c;
            if (Math.abs(a) < 1e-8) // Degenerate case
                return [];
            return [-b/a];
        }

        var D = b*b - 4*a*c;
        if (Math.abs(D) < 1e-8)
            return [-b/(2*a)];
        else if (D > 0)
            return [(-b+Math.sqrt(D))/(2*a), (-b-Math.sqrt(D))/(2*a)];
        return [];
    }

    // Convert to depressed cubic t^3+pt+q = 0 (subst x = t - b/3a)
    var p = (3*a*c - b*b)/(3*a*a);
    var q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/(27*a*a*a);
    var roots;

    if (Math.abs(p) < 1e-8) { // p = 0 -> t^3 = -q -> t = -q^1/3
        roots = [cuberoot(-q)];
    } else if (Math.abs(q) < 1e-8) { // q = 0 -> t^3 + pt = 0 -> t(t^2+p)=0
        roots = [0].concat(p < 0 ? [Math.sqrt(-p), -Math.sqrt(-p)] : []);
    } else {
        var D = q*q/4 + p*p*p/27;
        if (Math.abs(D) < 1e-8) {       // D = 0 -> two roots
            roots = [-1.5*q/p, 3*q/p];
        } else if (D > 0) {             // Only one real root
            var u = cuberoot(-q/2 - Math.sqrt(D));
            roots = [u - p/(3*u)];
        } else {                        // D < 0, three roots, but needs to use complex numbers/trigonometric solution
            var u = 2*Math.sqrt(-p/3);
            var t = Math.acos(3*q/p/u)/3;  // D < 0 implies p < 0 and acos argument in [-1..1]
            var k = 2*Math.PI/3;
            roots = [u*Math.cos(t), u*Math.cos(t-k), u*Math.cos(t-2*k)];
        }
    }

    // Convert back from depressed cubic
    for (var i = 0; i < roots.length; i++)
        roots[i] -= b/(3*a);

    return roots;
  }
  function cuberoot(x) {
    var y = Math.pow(Math.abs(x), 1/3);
    return x < 0 ? -y : y;
  }

  function stretch(dp, sp, p, min, max) { // integer
    var p1, p2, p3, p4;
    if (p < sp + 100) {
      p1 = min;
      p3 = min;
      p2 = sp + 100;
      p4 = dp + 100;
    } else {
      p1 = sp + 100;
      p3 = dp + 100;
      p2 = max;
      p4 = max;
    }
    return Math.floor(((p - p1) / (p2 - p1)) * (p4 - p3) + p3);
  }

  function getBoundingBox(strokes) { // minX, minY, maxX, maxY
    var a = new Object();
    a.minX = 200;
    a.minY = 200;
    a.maxX = 0;
    a.maxY = 0;
    for (var i = 0; i < strokes.length; i++) {
      if (strokes[i][0] == 0) { continue; }
      a.minX = Math.min(a.minX, strokes[i][3]);
      a.maxX = Math.max(a.maxX, strokes[i][3]);
      a.minY = Math.min(a.minY, strokes[i][4]);
      a.maxY = Math.max(a.maxY, strokes[i][4]);
      a.minX = Math.min(a.minX, strokes[i][5]);
      a.maxX = Math.max(a.maxX, strokes[i][5]);
      a.minY = Math.min(a.minY, strokes[i][6]);
      a.maxY = Math.max(a.maxY, strokes[i][6]);
      if (strokes[i][0] == 1) { continue; }
      if (strokes[i][0] == 99) { continue; }
      a.minX = Math.min(a.minX, strokes[i][7]);
      a.maxX = Math.max(a.maxX, strokes[i][7]);
      a.minY = Math.min(a.minY, strokes[i][8]);
      a.maxY = Math.max(a.maxY, strokes[i][8]);
      if (strokes[i][0] == 2) { continue; }
      if (strokes[i][0] == 3) { continue; }
      if (strokes[i][0] == 4) { continue; }
      a.minX = Math.min(a.minX, strokes[i][9]);
      a.maxX = Math.max(a.maxX, strokes[i][9]);
      a.minY = Math.min(a.minY, strokes[i][10]);
      a.maxY = Math.max(a.maxY, strokes[i][10]);
    }
    return a;
  }

  class PointMaker{
    constructor(x, y, dir, scale){
      this.x = x;
      this.y = y;
      this.dir = dir;
      this.scale = scale;
      if(!scale){
        this.scale = 1;
      }
      if(!dir){
      this.dir = DIR_POSX;//positive x
      }
    }
    setpos(x, y){
      this.x = x;
      this.y = y;
    }
    vec(x, y){ // void
      return [this.x + this.scale*this.dir.cos*x - this.scale*this.dir.sin*y,
              this.y + this.scale*this.dir.sin*x + this.scale*this.dir.cos*y]
    }
    setdir(dir){
      this.dir = dir;
    }
    setscale(scale){
      this.scale = scale;
    }
  }

  // ==ClosureCompiler==

  /**
   * Fit a Bezier curve to a (sub)set of digitized points.
   * Your code should not call this function directly. Use {@link fitCurve} instead.
   *
   * @param {Array<Array<Number>>} points - Array of digitized points, e.g. [[5,5],[5,50],[110,140],[210,160],[320,110]]
   * @param {Array<Number>} leftTangent - Unit tangent vector at start point
   * @param {Array<Number>} rightTangent - Unit tangent vector at end point
   * @param {Number} error - Tolerance, squared error between points and fitted curve
   * @returns {Array<Array<Array<Number>>>} Array of Bezier curves, where each element is [first-point, control-point-1, control-point-2, second-point] and points are [x, y]
   */

  function fitCubic_tang(points, tangents, error, progressCallback) {
      const MaxIterations = 20;   //Max times to try iterating (to find an acceptable curve)

      var bezCurve,               //Control points of fitted Bezier curve
          u,                      //Parameter values for point
          uPrime,                 //Improved parameter values
          maxError, prevErr,      //Maximum fitting error
          splitPoint, prevSplit,  //Point to split point set at if we need more than one curve
          //centerVector, toCenterTangent, fromCenterTangent,  //Unit tangent vector(s) at splitPoint
          beziers,                //Array of fitted Bezier curves if we need more than one curve
          dist, i;

      //console.log('fitCubic, ', points.length);
      const leftTangent = tangents[0];
      const rightTangent = maths.mulItems(tangents[points.length - 1], -1);
      //Use heuristic if region only has two points in it
      if (points.length === 2) {
          dist = maths.vectorLen(maths.subtract(points[0], points[1])) / 3.0;
          bezCurve = [
              points[0],
              maths.addArrays(points[0], maths.mulItems(leftTangent,  dist)),
              maths.addArrays(points[1], maths.mulItems(rightTangent, dist)),
              points[1]
          ];
          return [bezCurve];
      }


      //Parameterize points, and attempt to fit curve
      u = chordLengthParameterize(points);
      [bezCurve, maxError, splitPoint] = generateAndReport(points, u, u, leftTangent, rightTangent, progressCallback);

      if (maxError < error) {
          return [bezCurve];
      }
      //If error not too large, try some reparameterization and iteration
      if (maxError < (error*error)) {

          uPrime = u;
          prevErr = maxError;
          prevSplit = splitPoint;

          for (i = 0; i < MaxIterations; i++) {

              uPrime = reparameterize(bezCurve, points, uPrime);
              [bezCurve, maxError, splitPoint] = generateAndReport(points, u, uPrime, leftTangent, rightTangent, progressCallback);

              if (maxError < error) {
                  return [bezCurve];
              }
              //If the development of the fitted curve grinds to a halt,
              //we abort this attempt (and try a shorter curve):
              else if(splitPoint === prevSplit) {
                  let errChange = maxError/prevErr;
                  if((errChange > .9999) && (errChange < 1.0001)) {
                      break;
                  }
              }

              prevErr = maxError;
              prevSplit = splitPoint;
          }
      }

      //Fitting failed -- split at max error point and fit recursively
      beziers = [];

      //To and from need to point in opposite directions:
      /*
      Note: An alternative to this "divide and conquer" recursion could be to always
            let new curve segments start by trying to go all the way to the end,
            instead of only to the end of the current subdivided polyline.
            That might let many segments fit a few points more, reducing the number of total segments.

            However, a few tests have shown that the segment reduction is insignificant
            (240 pts, 100 err: 25 curves vs 27 curves. 140 pts, 100 err: 17 curves on both),
            and the results take twice as many steps and milliseconds to finish,
            without looking any better than what we already have.
      */
      beziers = beziers.concat(fitCubic_tang(points.slice(0, splitPoint + 1), tangents.slice(0, splitPoint + 1), error, progressCallback));
      beziers = beziers.concat(fitCubic_tang(points.slice(splitPoint),        tangents.slice(splitPoint)       , error, progressCallback));
      return beziers;
  }
  function generateAndReport(points, paramsOrig, paramsPrime, leftTangent, rightTangent, progressCallback) {
      var bezCurve, maxError, splitPoint;

      bezCurve = generateBezier(points, paramsPrime, leftTangent, rightTangent);
      //Find max deviation of points to fitted curve.
      //Here we always use the original parameters (from chordLengthParameterize()),
      //because we need to compare the current curve to the actual source polyline,
      //and not the currently iterated parameters which reparameterize() & generateBezier() use,
      //as those have probably drifted far away and may no longer be in ascending order.
      [maxError, splitPoint] = computeMaxError(points, bezCurve, paramsOrig);

      if(progressCallback) {
          progressCallback({
              bez: bezCurve,
              points: points,
              params: paramsOrig,
              maxErr: maxError,
              maxPoint: splitPoint,
          });
      }

      return [bezCurve, maxError, splitPoint];
  }

  /**
   * Use least-squares method to find Bezier control points for region.
   *
   * @param {Array<Array<Number>>} points - Array of digitized points
   * @param {Array<Number>} parameters - Parameter values for region
   * @param {Array<Number>} leftTangent - Unit tangent vector at start point
   * @param {Array<Number>} rightTangent - Unit tangent vector at end point
   * @returns {Array<Array<Number>>} Approximated Bezier curve: [first-point, control-point-1, control-point-2, second-point] where points are [x, y]
   */
  function generateBezier(points, parameters, leftTangent, rightTangent) {
      var bezCurve,                       //Bezier curve ctl pts
          A, a,                           //Precomputed rhs for eqn
          C, X,                           //Matrices C & X
          det_C0_C1, det_C0_X, det_X_C1,  //Determinants of matrices
          alpha_l, alpha_r,               //Alpha values, left and right

          epsilon, segLength,
          i, len, tmp, u, ux,
          firstPoint = points[0],
          lastPoint = points[points.length-1];

      bezCurve = [firstPoint, null, null, lastPoint];
      //console.log('gb', parameters.length);

      //Compute the A's
      A = maths.zeros_Xx2x2(parameters.length);
      for (i = 0, len = parameters.length; i < len; i++) {
          u = parameters[i];
          ux = 1 - u;
          a = A[i];

          a[0] = maths.mulItems(leftTangent,  3 * u  * (ux*ux));
          a[1] = maths.mulItems(rightTangent, 3 * ux * (u*u));
      }

      //Create the C and X matrices
      C = [[0,0], [0,0]];
      X = [0,0];
      for (i = 0, len = points.length; i < len; i++) {
          u = parameters[i];
          a = A[i];

          C[0][0] += maths.dot(a[0], a[0]);
          C[0][1] += maths.dot(a[0], a[1]);
          C[1][0] += maths.dot(a[0], a[1]);
          C[1][1] += maths.dot(a[1], a[1]);

          tmp = maths.subtract(points[i], bezier.q([firstPoint, firstPoint, lastPoint, lastPoint], u));

          X[0] += maths.dot(a[0], tmp);
          X[1] += maths.dot(a[1], tmp);
      }

      //Compute the determinants of C and X
      det_C0_C1 = (C[0][0] * C[1][1]) - (C[1][0] * C[0][1]);
      det_C0_X  = (C[0][0] * X[1]   ) - (C[1][0] * X[0]   );
      det_X_C1  = (X[0]    * C[1][1]) - (X[1]    * C[0][1]);

      //Finally, derive alpha values
      alpha_l = det_C0_C1 === 0 ? 0 : det_X_C1 / det_C0_C1;
      alpha_r = det_C0_C1 === 0 ? 0 : det_C0_X / det_C0_C1;

      //If alpha negative, use the Wu/Barsky heuristic (see text).
      //If alpha is 0, you get coincident control points that lead to
      //divide by zero in any subsequent NewtonRaphsonRootFind() call.
      segLength = maths.vectorLen(maths.subtract(firstPoint, lastPoint));
      epsilon = 1.0e-6 * segLength;
      if (alpha_l < epsilon || alpha_r < epsilon) {
          //Fall back on standard (probably inaccurate) formula, and subdivide further if needed.
          bezCurve[1] = maths.addArrays(firstPoint, maths.mulItems(leftTangent,  segLength / 3.0));
          bezCurve[2] = maths.addArrays(lastPoint,  maths.mulItems(rightTangent, segLength / 3.0));
      } else {
          //First and last control points of the Bezier curve are
          //positioned exactly at the first and last data points
          //Control points 1 and 2 are positioned an alpha distance out
          //on the tangent vectors, left and right, respectively
          bezCurve[1] = maths.addArrays(firstPoint, maths.mulItems(leftTangent,  alpha_l));
          bezCurve[2] = maths.addArrays(lastPoint,  maths.mulItems(rightTangent, alpha_r));
      }

      return bezCurve;
  }
  /**
   * Given set of points and their parameterization, try to find a better parameterization.
   *
   * @param {Array<Array<Number>>} bezier - Current fitted curve
   * @param {Array<Array<Number>>} points - Array of digitized points
   * @param {Array<Number>} parameters - Current parameter values
   * @returns {Array<Number>} New parameter values
   */
  function reparameterize(bezier, points, parameters) {
      /*
      var j, len, point, results, u;
      results = [];
      for (j = 0, len = points.length; j < len; j++) {
          point = points[j], u = parameters[j];

          results.push(newtonRaphsonRootFind(bezier, point, u));
      }
      return results;
      //*/
      return parameters.map((p, i) => newtonRaphsonRootFind(bezier, points[i], p));
  }
  /**
   * Use Newton-Raphson iteration to find better root.
   *
   * @param {Array<Array<Number>>} bez - Current fitted curve
   * @param {Array<Number>} point - Digitized point
   * @param {Number} u - Parameter value for "P"
   * @returns {Number} New u
   */
  function newtonRaphsonRootFind(bez, point, u) {
      /*
          Newton's root finding algorithm calculates f(x)=0 by reiterating
          x_n+1 = x_n - f(x_n)/f'(x_n)
          We are trying to find curve parameter u for some point p that minimizes
          the distance from that point to the curve. Distance point to curve is d=q(u)-p.
          At minimum distance the point is perpendicular to the curve.
          We are solving
          f = q(u)-p * q'(u) = 0
          with
          f' = q'(u) * q'(u) + q(u)-p * q''(u)
          gives
          u_n+1 = u_n - |q(u_n)-p * q'(u_n)| / |q'(u_n)**2 + q(u_n)-p * q''(u_n)|
      */

      var d = maths.subtract(bezier.q(bez, u), point),
          qprime = bezier.qprime(bez, u),
          numerator = maths.mulMatrix(d, qprime),
          denominator = maths.sum(maths.squareItems(qprime)) + 2 * maths.mulMatrix(d, bezier.qprimeprime(bez, u));

      if (denominator === 0) {
          return u;
      } else {
          return u - (numerator/denominator);
      }
  }
  /**
   * Assign parameter values to digitized points using relative distances between points.
   *
   * @param {Array<Array<Number>>} points - Array of digitized points
   * @returns {Array<Number>} Parameter values
   */
  function chordLengthParameterize(points) {
      var u = [], currU, prevU, prevP;

      points.forEach((p, i) => {
          currU = i ? prevU + maths.vectorLen(maths.subtract(p, prevP))
                    : 0;
          u.push(currU);

          prevU = currU;
          prevP = p;
      });
      u = u.map(x => x/prevU);

      return u;
  }
  /**
   * Find the maximum squared distance of digitized points to fitted curve.
   *
   * @param {Array<Array<Number>>} points - Array of digitized points
   * @param {Array<Array<Number>>} bez - Fitted curve
   * @param {Array<Number>} parameters - Parameterization of points
   * @returns {Array<Number>} Maximum error (squared) and point of max error
   */
  function computeMaxError(points, bez, parameters) {
      var dist,       //Current error
          maxDist,    //Maximum error
          splitPoint, //Point of maximum error
          v,          //Vector from point to curve
          i, count, point, t;

      maxDist = 0;
      splitPoint = points.length / 2;

      const t_distMap = mapTtoRelativeDistances(bez, 10);

      for (i = 0, count = points.length; i < count; i++) {
          point = points[i];
          //Find 't' for a point on the bez curve that's as close to 'point' as possible:
          t = find_t(bez, parameters[i], t_distMap, 10);

          v = maths.subtract(bezier.q(bez, t), point);
          dist = v[0]*v[0] + v[1]*v[1];

          if (dist > maxDist) {
              maxDist = dist;
              splitPoint = i;
          }
      }

      return [maxDist, splitPoint];
  }
  //Sample 't's and map them to relative distances along the curve:
  var mapTtoRelativeDistances = function (bez, B_parts) {
      var B_t_curr;
      var B_t_dist = [0];
      var B_t_prev = bez[0];
      var sumLen = 0;

      for (var i=1; i<=B_parts; i++) {
        B_t_curr = bezier.q(bez, i/B_parts);

        sumLen += maths.vectorLen(maths.subtract(B_t_curr, B_t_prev));

        B_t_dist.push(sumLen);
        B_t_prev = B_t_curr;
      }

      //Normalize B_length to the same interval as the parameter distances; 0 to 1:
      B_t_dist = B_t_dist.map(x => x/sumLen);
      return B_t_dist;
  };

  function find_t(bez, param, t_distMap, B_parts) {
      if(param < 0) { return 0; }
      if(param > 1) { return 1; }

      /*
          'param' is a value between 0 and 1 telling us the relative position
          of a point on the source polyline (linearly from the start (0) to the end (1)).
          To see if a given curve - 'bez' - is a close approximation of the polyline,
          we compare such a poly-point to the point on the curve that's the same
          relative distance along the curve's length.

          But finding that curve-point takes a little work:
          There is a function "B(t)" to find points along a curve from the parametric parameter 't'
          (also relative from 0 to 1: http://stackoverflow.com/a/32841764/1869660
                                      http://pomax.github.io/bezierinfo/#explanation),
          but 't' isn't linear by length (http://gamedev.stackexchange.com/questions/105230).

          So, we sample some points along the curve using a handful of values for 't'.
          Then, we calculate the length between those samples via plain euclidean distance;
          B(t) concentrates the points around sharp turns, so this should give us a good-enough outline of the curve.
          Thus, for a given relative distance ('param'), we can now find an upper and lower value
          for the corresponding 't' by searching through those sampled distances.
          Finally, we just use linear interpolation to find a better value for the exact 't'.

          More info:
              http://gamedev.stackexchange.com/questions/105230/points-evenly-spaced-along-a-bezier-curve
              http://stackoverflow.com/questions/29438398/cheap-way-of-calculating-cubic-bezier-length
              http://steve.hollasch.net/cgindex/curves/cbezarclen.html
              https://github.com/retuxx/tinyspline
      */
      var lenMax, lenMin, tMax, tMin, t;

      //Find the two t-s that the current param distance lies between,
      //and then interpolate a somewhat accurate value for the exact t:
      for(var i = 1; i <= B_parts; i++) {

          if(param <= t_distMap[i]) {
              tMin   = (i-1) / B_parts;
              tMax   = i / B_parts;
              lenMin = t_distMap[i-1];
              lenMax = t_distMap[i];

              t = (param-lenMin)/(lenMax-lenMin) * (tMax-tMin) + tMin;
              break;
          }
      }
      return t;
  }

  /*
      Simplified versions of what we need from math.js
      Optimized for our input, which is only numbers and 1x2 arrays (i.e. [x, y] coordinates).
  */
  class maths {
      //zeros = logAndRun(math.zeros);
      static zeros_Xx2x2(x) {
          var zs = [];
          while(x--) { zs.push([0,0]); }
          return zs
      }

      //multiply = logAndRun(math.multiply);
      static mulItems(items, multiplier) {
          return items.map(x => x*multiplier);
      }
      static mulMatrix(m1, m2) {
          //https://en.wikipedia.org/wiki/Matrix_multiplication#Matrix_product_.28two_matrices.29
          //Simplified to only handle 1-dimensional matrices (i.e. arrays) of equal length:
           return m1.reduce((sum,x1,i) => sum + (x1*m2[i]), 0);
      }

      //Only used to subract to points (or at least arrays):
      //  subtract = logAndRun(math.subtract);
      static subtract(arr1, arr2) {
          return arr1.map((x1, i) => x1 - arr2[i]);
      }

      //add = logAndRun(math.add);
      static addArrays(arr1, arr2) {
          return arr1.map((x1, i) => x1 + arr2[i]);
      }
      static addItems(items, addition) {
          return items.map(x => x+addition);
      }

      //var sum = logAndRun(math.sum);
      static sum(items) {
          return items.reduce((sum,x) => sum + x);
      }

      //chain = math.chain;

      //Only used on two arrays. The dot product is equal to the matrix product in this case:
      //  dot = logAndRun(math.dot);
      static dot(m1, m2) {
          return maths.mulMatrix(m1, m2);
      }

      //https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm
      //  var norm = logAndRun(math.norm);
      static vectorLen(v) {
          return Math.hypot(...v);
      }

      //math.divide = logAndRun(math.divide);
      static divItems(items, divisor) {
          return items.map(x => x/divisor);
      }

      //var dotPow = logAndRun(math.dotPow);
      static squareItems(items) {
          return items.map(x => x*x);
      }

      static normalize(v) {
          return this.divItems(v, this.vectorLen(v));
      }

      //Math.pow = logAndRun(Math.pow);
  }


  class bezier {
      //Evaluates cubic bezier at t, return point
      static q(ctrlPoly, t) {
          var tx = 1.0 - t;
          var pA = maths.mulItems( ctrlPoly[0],      tx * tx * tx ),
              pB = maths.mulItems( ctrlPoly[1],  3 * tx * tx *  t ),
              pC = maths.mulItems( ctrlPoly[2],  3 * tx *  t *  t ),
              pD = maths.mulItems( ctrlPoly[3],       t *  t *  t );
          return maths.addArrays(maths.addArrays(pA, pB), maths.addArrays(pC, pD));
      }

      //Evaluates cubic bezier first derivative at t, return point
      static qprime(ctrlPoly, t) {
          var tx = 1.0 - t;
          var pA = maths.mulItems( maths.subtract(ctrlPoly[1], ctrlPoly[0]),  3 * tx * tx ),
              pB = maths.mulItems( maths.subtract(ctrlPoly[2], ctrlPoly[1]),  6 * tx *  t ),
              pC = maths.mulItems( maths.subtract(ctrlPoly[3], ctrlPoly[2]),  3 *  t *  t );
          return maths.addArrays(maths.addArrays(pA, pB), pC);
      }

      //Evaluates cubic bezier second derivative at t, return point
      static qprimeprime(ctrlPoly, t) {
          return maths.addArrays(maths.mulItems( maths.addArrays(maths.subtract(ctrlPoly[2], maths.mulItems(ctrlPoly[1], 2)), ctrlPoly[0]),  6 * (1.0 - t) ),
                                 maths.mulItems( maths.addArrays(maths.subtract(ctrlPoly[3], maths.mulItems(ctrlPoly[2], 2)), ctrlPoly[1]),  6 *        t  ));
      }
  }

  //module.exports = fitCurve;

  const bezier_steps = 200;
  class Bezier{
    static generalBezier(x_fun, y_fun, dx_fun, dy_fun, width_func, width_func_d){
      var a1 = [];
      var a2 = [];
      var tang1 = [];
      var tang2 = [];
      for (var tt = 0; tt <= bezier_steps; tt++) {
        const t = tt / bezier_steps;
        const x = x_fun(t);
        const y = y_fun(t);
        const vx = dx_fun(t);
        const vy = dy_fun(t);

        let [ia, ib] = unit_normal_vector(vx, vy);
        const deltad = width_func(t);
        ia = ia * deltad;
        ib = ib * deltad;
        
        const rad = get_rad(vx, vy);
        const velocity = Math.sqrt(vx*vx+vy*vy);
        const width_rad = Math.atan(width_func_d(t)/velocity);
        a1.push([x - ia, y - ib]);
        a2.push([x + ia, y + ib]);
        tang1.push(rad_to_vector(rad-width_rad));
        tang2.push(rad_to_vector(rad+width_rad-Math.PI));
      }
      const bez1 = fitCubic_tang(a1, tang1, 0.03);
      const bez2 = fitCubic_tang(a2.reverse(), tang2.reverse(), 0.03);
      //const bez1 = fitCurve(a1, 0.01);
      //const bez2 = fitCurve(a2.reverse(), 0.01);
      return [bez1, bez2];
    }
    static generalBezier2(x_fun, y_fun, dx_fun, dy_fun, width_func, width_func_d, dir_func){
      //offset vector (ia, ib) is calculated with dir_func
      var a1 = [];
      var a2 = [];
      var tang1 = [];
      var tang2 = [];
      for (var tt = 0; tt <= bezier_steps; tt++) {
        const t = tt / bezier_steps;
        const x = x_fun(t);
        const y = y_fun(t);
        const vx = dx_fun(t);
        const vy = dy_fun(t);

        let [ia, ib] = dir_func(t);
        const deltad = width_func(t);
        ia = ia * deltad;
        ib = ib * deltad;
        
        const rad = get_rad(vx, vy);
        const velocity = Math.sqrt(vx*vx+vy*vy);
        const width_rad = Math.atan(width_func_d(t)/velocity);
        a1.push([x - ia, y - ib]);
        a2.push([x + ia, y + ib]);
        tang1.push(rad_to_vector(rad-width_rad));
        tang2.push(rad_to_vector(rad+width_rad-Math.PI));
      }
      const bez1 = fitCubic_tang(a1, tang1, 0.03);
      const bez2 = fitCubic_tang(a2.reverse(), tang2.reverse(), 0.03);
      //const bez1 = fitCurve(a1, 0.03);
      //const bez2 = fitCurve(a2.reverse(), 0.03);
      return [bez1, bez2];
    }

    static qBezier(x1, y1, sx, sy, x2, y2, width_func, width_func_d){
      const x_fun = t => ((1.0 - t) * (1.0 - t) * x1 + 2.0 * t * (1.0 - t) * sx + t * t * x2);
      const y_fun = t => ((1.0 - t) * (1.0 - t) * y1 + 2.0 * t * (1.0 - t) * sy + t * t * y2);
      const dx_fun = t => (x1 - 2.0 * sx + x2) * 2.0 * t + (-2.0 * x1 + 2.0 * sx);
      const dy_fun = t => (y1 - 2.0 * sy + y2) * 2.0 * t + (-2.0 * y1 + 2.0 * sy);
      return this.generalBezier(x_fun, y_fun, dx_fun, dy_fun, width_func, width_func_d);
    }
    static qBezier2(x1, y1, sx, sy, x2, y2, width_func, width_func_d){
      //similar to qBezier(), but the direction changes at a constant speed (not decided by normal vector)
      const x_fun = t => ((1.0 - t) * (1.0 - t) * x1 + 2.0 * t * (1.0 - t) * sx + t * t * x2);
      const y_fun = t => ((1.0 - t) * (1.0 - t) * y1 + 2.0 * t * (1.0 - t) * sy + t * t * y2);
      const dx_fun = t => (x1 - 2.0 * sx + x2) * 2.0 * t + (-2.0 * x1 + 2.0 * sx);
      const dy_fun = t => (y1 - 2.0 * sy + y2) * 2.0 * t + (-2.0 * y1 + 2.0 * sy);
      const rad_begin = Math.atan2(sy-y1, sx-x1);
      var rad_end = Math.atan2(y2-sy, x2-sx);
      if(rad_end - rad_begin > Math.PI) {
        rad_end -= Math.PI*2;
      }else if(rad_begin - rad_end > Math.PI){
        rad_end += Math.PI*2;
      }
      const dir_func = t => rad_to_vector((1-t)*rad_begin+t*rad_end+Math.PI/2);
      return this.generalBezier2(x_fun, y_fun, dx_fun, dy_fun, width_func, width_func_d, dir_func);
    }
    
    static cBezier(x1, y1, sx1, sy1, sx2, sy2, x2, y2, width_func, width_func_d){
      const x_fun = t => (1.0 - t) * (1.0 - t) * (1.0 - t) * x1 + 3.0 * t * (1.0 - t) * (1.0 - t) * sx1 + 3 * t * t * (1.0 - t) * sx2 + t * t * t * x2;
      const y_fun = t => (1.0 - t) * (1.0 - t) * (1.0 - t) * y1 + 3.0 * t * (1.0 - t) * (1.0 - t) * sy1 + 3 * t * t * (1.0 - t) * sy2 + t * t * t * y2;
      const dx_fun = t => t * t * (-3 * x1 + 9 * sx1 + -9 * sx2 + 3 * x2) + t * (6 * x1 + -12 * sx1 + 6 * sx2) + -3 * x1 + 3 * sx1;
      const dy_fun = t => t * t * (-3 * y1 + 9 * sy1 + -9 * sy2 + 3 * y2) + t * (6 * y1 + -12 * sy1 + 6 * sy2) + -3 * y1 + 3 * sy1;
      return this.generalBezier(x_fun, y_fun, dx_fun, dy_fun, width_func, width_func_d);
    }
    
    static slantBezier(x_fun, y_fun, dx_fun, dy_fun, width_func, width_func_d, dir_x, dir_y){
      var a1 = [];
      var a2 = [];
      var tang1 = [];
      var tang2 = [];
      let [ia, ib] = unit_normal_vector(dir_x, dir_y);
      let len = Math.sqrt(dir_x*dir_x+dir_y*dir_y);
      let ex = dir_x/len;
      let ey = dir_y/len;
      
      for (var tt = 0; tt <= bezier_steps; tt++) {
        const t = tt / bezier_steps;
        const x = x_fun(t);
        const y = y_fun(t);
        const vx = dx_fun(t);
        const vy = dy_fun(t);

        const deltad = width_func(t);
        
        const velocity = (dir_x*vx+dir_y*vy)/Math.sqrt(dir_x*dir_x+dir_y*dir_y);
        const width_tan = width_func_d(t)/velocity;
        const bez_tan = (dir_x*vy-dir_y*vx)/(dir_x*vx+dir_y*vy);
        const rad1 = Math.atan(bez_tan-width_tan);
        const rad2 = Math.atan(bez_tan+width_tan);
        a1.push([x - ia * deltad, y - ib * deltad]);
        a2.push([x + ia * deltad, y + ib * deltad]);
        tang1.push([ Math.cos(rad1)*ex-Math.sin(rad1)*ey,  Math.sin(rad1)*ex+Math.cos(rad1)*ey]);
        tang2.push([-Math.cos(rad2)*ex+Math.sin(rad2)*ey, -Math.sin(rad2)*ex-Math.cos(rad2)*ey]);
      }
      const bez1 = fitCubic_tang(a1, tang1, 0.03);
      const bez2 = fitCubic_tang(a2.reverse(), tang2.reverse(), 0.03);
      //const bez1 = fitCurve(a1, 0.01);
      //const bez2 = fitCurve(a2.reverse(), 0.01);
      return [bez1, bez2];
    }

    static qBezier_slant(x1, y1, sx, sy, x2, y2, width_func, width_func_d){
      const x_fun = t => ((1.0 - t) * (1.0 - t) * x1 + 2.0 * t * (1.0 - t) * sx + t * t * x2);
      const y_fun = t => ((1.0 - t) * (1.0 - t) * y1 + 2.0 * t * (1.0 - t) * sy + t * t * y2);
      const dir_x = x2 - x1;
      const dir_y = y2 - y1;
      const dx_fun = t => (x1 - 2.0 * sx + x2) * 2.0 * t + (-2.0 * x1 + 2.0 * sx);
      const dy_fun = t => (y1 - 2.0 * sy + y2) * 2.0 * t + (-2.0 * y1 + 2.0 * sy);
      return this.slantBezier(x_fun, y_fun, dx_fun, dy_fun, width_func, width_func_d, dir_x, dir_y);
    }
    static cBezier_slant(x1, y1, sx1, sy1, sx2, sy2, x2, y2, width_func, width_func_d){
      const x_fun = t => (1.0 - t) * (1.0 - t) * (1.0 - t) * x1 + 3.0 * t * (1.0 - t) * (1.0 - t) * sx1 + 3 * t * t * (1.0 - t) * sx2 + t * t * t * x2;
      const y_fun = t => (1.0 - t) * (1.0 - t) * (1.0 - t) * y1 + 3.0 * t * (1.0 - t) * (1.0 - t) * sy1 + 3 * t * t * (1.0 - t) * sy2 + t * t * t * y2;
      const dx_fun = t => t * t * (-3 * x1 + 9 * sx1 + -9 * sx2 + 3 * x2) + t * (6 * x1 + -12 * sx1 + 6 * sx2) + -3 * x1 + 3 * sx1;
      const dy_fun = t => t * t * (-3 * y1 + 9 * sy1 + -9 * sy2 + 3 * y2) + t * (6 * y1 + -12 * sy1 + 6 * sy2) + -3 * y1 + 3 * sy1;
      const dir_x = x2 - x1;
      const dir_y = y2 - y1;
      return this.slantBezier(x_fun, y_fun, dx_fun, dy_fun, width_func, width_func_d, dir_x, dir_y);
    }
    
   static bez_to_poly(bez){
      var poly = new Polygon();
      poly.push(bez[0][0][0], bez[0][0][1]);
      for(let bez1 of bez){
        poly.push(bez1[1][0], bez1[1][1], 2);
        poly.push(bez1[2][0], bez1[2][1], 2);
        poly.push(bez1[3][0], bez1[3][1]);
      }
      return poly;
    }
  }

  class FontCanvas {
    constructor() {
      this.polygons = new Polygons();
    }
    getPolygons() {
      return this.polygons;
    }
    addPolygon(poly) {
      this.polygons.push(poly);
    }

    //flipping and rotating
    //polygons in the rectangle (x1, y1, x2, y2) are flipped or rotated
    flip_left_right(x1, y1, x2, y2) {
      for (var i = 0; i < this.polygons.array.length; i++) {
        var inside = true;
        for (var j = 0; j < this.polygons.array[i].array.length; j++) {
          if (x1 > this.polygons.array[i].array[j].x || this.polygons.array[i].array[j].x > x2 ||
            y1 > this.polygons.array[i].array[j].y || this.polygons.array[i].array[j].y > y2) {
            inside = false;
          }
        }
        if (inside) {
          for (var j = 0; j < this.polygons.array[i].array.length; j++) {
            this.polygons.array[i].array[j].x = x2 - (this.polygons.array[i].array[j].x - x1);
            this.polygons.array[i].array[j].x = Math.floor(this.polygons.array[i].array[j].x * 10) / 10;
          }
        }
      }
    }

    flip_up_down(x1, y1, x2, y2) {
      for (var i = 0; i < this.polygons.array.length; i++) {
        var inside = true;
        for (var j = 0; j < this.polygons.array[i].array.length; j++) {
          if (x1 > this.polygons.array[i].array[j].x || this.polygons.array[i].array[j].x > x2 ||
            y1 > this.polygons.array[i].array[j].y || this.polygons.array[i].array[j].y > y2) {
            inside = false;
          }
        }
        if (inside) {
          for (var j = 0; j < this.polygons.array[i].array.length; j++) {
            this.polygons.array[i].array[j].y = y2 - (this.polygons.array[i].array[j].y - y1);
            this.polygons.array[i].array[j].y = Math.floor(this.polygons.array[i].array[j].y * 10) / 10;
          }
        }
      }
    }
    
    rotate90(x1, y1, x2, y2) {
      for (var i = 0; i < this.polygons.array.length; i++) {
        var inside = true;
        for (var j = 0; j < this.polygons.array[i].array.length; j++) {
          if (x1 > this.polygons.array[i].array[j].x || this.polygons.array[i].array[j].x > x2 ||
            y1 > this.polygons.array[i].array[j].y || this.polygons.array[i].array[j].y > y2) {
            inside = false;
          }
        }
        if (inside) {
          for (var j = 0; j < this.polygons.array[i].array.length; j++) {
            var x = this.polygons.array[i].array[j].x;
            var y = this.polygons.array[i].array[j].y;
            this.polygons.array[i].array[j].x = x1 + (y2 - y);
            this.polygons.array[i].array[j].y = y1 + (x - x1);
            this.polygons.array[i].array[j].x = Math.floor(this.polygons.array[i].array[j].x * 10) / 10;
            this.polygons.array[i].array[j].y = Math.floor(this.polygons.array[i].array[j].y * 10) / 10;
          }
        }
      }
    }

    rotate180(x1, y1, x2, y2) {
      for (var i = 0; i < this.polygons.array.length; i++) {
        var inside = true;
        for (var j = 0; j < this.polygons.array[i].array.length; j++) {
          if (x1 > this.polygons.array[i].array[j].x || this.polygons.array[i].array[j].x > x2 ||
            y1 > this.polygons.array[i].array[j].y || this.polygons.array[i].array[j].y > y2) {
            inside = false;
          }
        }
        if (inside) {
          for (var j = 0; j < this.polygons.array[i].array.length; j++) {
            var x = this.polygons.array[i].array[j].x;
            var y = this.polygons.array[i].array[j].y;
            this.polygons.array[i].array[j].x = x2 - (x - x1);
            this.polygons.array[i].array[j].y = y2 - (y - y1);
            this.polygons.array[i].array[j].x = Math.floor(this.polygons.array[i].array[j].x * 10) / 10;
            this.polygons.array[i].array[j].y = Math.floor(this.polygons.array[i].array[j].y * 10) / 10;
          }
        }
      }
    }

    rotate270(x1, y1, x2, y2) {
      for (var i = 0; i < this.polygons.array.length; i++) {
        var inside = true;
        for (var j = 0; j < this.polygons.array[i].array.length; j++) {
          if (x1 > this.polygons.array[i].array[j].x || this.polygons.array[i].array[j].x > x2 ||
            y1 > this.polygons.array[i].array[j].y || this.polygons.array[i].array[j].y > y2) {
            inside = false;
          }
        }
        if (inside) {
          for (var j = 0; j < this.polygons.array[i].array.length; j++) {
            var x = this.polygons.array[i].array[j].x;
            var y = this.polygons.array[i].array[j].y;
            this.polygons.array[i].array[j].x = x1 + (y - y1);
            this.polygons.array[i].array[j].y = y2 - (x - x1);
            this.polygons.array[i].array[j].x = Math.floor(this.polygons.array[i].array[j].x * 10) / 10;
            this.polygons.array[i].array[j].y = Math.floor(this.polygons.array[i].array[j].y * 10) / 10;
          }
        }
      }
    }

    drawUpperLeftCorner(x1, y1, dir, kMinWidthT) {
      var p = new PointMaker(x1, y1, dir, kMinWidthT);
      var poly = new Polygon();
      poly.push2(p.vec(0, 1));
      poly.push2(p.vec(0, -1));
      poly.push2(p.vec(-1, 1));
      poly.reverse();
      this.polygons.push(poly);
    }

    drawUpperRightCorner(x1, y1, kMinWidthT, kagekMinWidthY, kagekWidth) {
      var p = new PointMaker(x1, y1);
      var poly = new Polygon();
      poly.push2(p.vec(-kMinWidthT,-kagekMinWidthY));
      poly.push2(p.vec(0,-kagekMinWidthY - kagekWidth));
      poly.push2(p.vec(kMinWidthT+kagekWidth,kagekMinWidthY));
      poly.push2(p.vec(kMinWidthT,kMinWidthT*0.8), 2);
      poly.push2(p.vec(0,kMinWidthT*1.2), 2);
      poly.push2(p.vec(-kMinWidthT,kMinWidthT*1.3));
      this.polygons.push(poly);
    }
    
    drawUpperRightCorner_straight_v(x1, y1, kMinWidthT, kagekMinWidthY, kagekWidth) {//vertical straight line
      var p = new PointMaker(x1, y1);
      var poly = new Polygon();
      poly.push2(p.vec(-kMinWidthT,-kagekMinWidthY));
      poly.push2(p.vec(0,-kagekMinWidthY - kagekWidth));
      poly.push2(p.vec(kMinWidthT+kagekWidth,kagekMinWidthY));
      poly.push2(p.vec(kMinWidthT,kMinWidthT));
      poly.push2(p.vec(-kMinWidthT,0));
      this.polygons.push(poly);
    }

    drawOpenBegin_straight(x1, y1, kMinWidthT, kagekMinWidthY, dir) {
      const rad_offset = Math.atan(kagekMinWidthY*0.5/kMinWidthT);
      //const rad_offset = (0.1*Math.PI);
      var poly = new Polygon();
      let p1 = new PointMaker(x1, y1, dir);
      let[x, y] = p1.vec(kagekMinWidthY*0.5, -kMinWidthT);
      const offs_sin = Math.sin(rad_offset);
      const offs_cos = Math.cos(rad_offset);
      //let[x, y] = p1.vec(kMinWidthT*offs_sin/offs_cos, -kMinWidthT);
      const new_dir = {sin: dir.sin*offs_cos+dir.cos*offs_sin, cos: dir.cos*offs_cos-dir.sin*offs_sin};
      let p2 = new PointMaker(x, y, new_dir, kMinWidthT);
      poly.push2(p2.vec(0, 0));
      poly.push2(p2.vec(0, -1.4), 2);
      poly.push2(p2.vec(0.6, -1.4), 2);
      poly.push2(p2.vec(2.0, 1.0));
      this.polygons.push(poly);
    }

    drawOpenBegin_curve_down2(x1, y1, dir, kMinWidthT, rad_offset){
      var poly = new Polygon();
      let p1 = new PointMaker(x1, y1, dir);
      poly.push2(p1.vec(0, kMinWidthT));

      let p2 = new PointMaker();
      const offs_sin = Math.sin(rad_offset);
      const offs_cos = Math.cos(rad_offset);
      p2.setscale(kMinWidthT);
      p2.setdir({sin: dir.sin*offs_cos+dir.cos*offs_sin, cos: dir.cos*offs_cos-dir.sin*offs_sin});
      if(rad_offset>0){
        poly.push2(p1.vec(-(offs_sin/offs_cos)*2*kMinWidthT, kMinWidthT));
        let[x, y] = p1.vec(0, -kMinWidthT);
        poly.push(x, y);
        p2.setpos(x, y);
        let[x2, y2] = p2.vec(0, -rad_offset);
        p2.setpos(x2, y2);
      }else {
        let[x, y] = p1.vec(0+2*kMinWidthT*offs_sin, kMinWidthT-2*kMinWidthT*offs_cos);
        p2.setpos(x, y);
      }

      poly.push2(p2.vec(0, 0));
      poly.push2(p2.vec(0, -0.8), 2);
      poly.push2(p2.vec(0.6, -0.8), 2);
      poly.push2(p2.vec(1.8, 1.0));
      this.polygons.push(poly);
    }

    drawOpenBegin_curve_down(x1, y1, dir, kMinWidthT, kagekMinWidthY){
      const rad = Math.atan2(Math.abs(dir.sin), Math.abs(dir.cos));
      var rad_offset;
      if(rad > Math.PI * 0.2){//36 degrees
        rad_offset = (0.1*Math.PI)*(rad-Math.PI * 0.2)/(Math.PI*0.3);
      }else {
        rad_offset = (-0.25*Math.PI)*(Math.PI*0.2-rad)/(Math.PI*0.2);
      }
      this.drawOpenBegin_curve_down2(x1, y1, dir, kMinWidthT, rad_offset);
    }

    drawOpenBegin_curve_up(x1, y1, dir, kMinWidthT, kagekMinWidthY) {
      var poly = new Polygon();
      let p1 = new PointMaker(x1, y1, dir);
      const offs_sin = Math.sin(-0.1*Math.PI);
      const offs_cos = Math.cos(-0.1*Math.PI);
      poly.push2(p1.vec(0, -kMinWidthT));
      poly.push2(p1.vec((offs_sin/offs_cos)*2*kMinWidthT, -kMinWidthT));

      let p2 = new PointMaker();
      p2.setscale(kMinWidthT);
      p2.setdir({sin: dir.sin*offs_cos+dir.cos*offs_sin, cos: dir.cos*offs_cos-dir.sin*offs_sin});
      let[x, y] = p1.vec(0, kMinWidthT);
      p2.setpos(x, y);

      poly.push2(p2.vec(0, 0));
      poly.push2(p2.vec(0, 0.8), 2);
      poly.push2(p2.vec(0.6, 0.8), 2);
      poly.push2(p2.vec(1.8, -1.0));
      poly.reverse();
      this.polygons.push(poly);
    }

    drawLowerRightHT_v(x2, y2, kMinWidthT, kagekMinWidthY) {//for T design
      var poly = new Polygon();
      poly.push(x2, y2 + kagekMinWidthY);
      poly.push(x2 + kMinWidthT, y2 - kagekMinWidthY * 3);
      poly.push(x2 + kMinWidthT * 2, y2 - kagekMinWidthY);
      poly.push(x2 + kMinWidthT * 2, y2 + kagekMinWidthY);
      this.polygons.push(poly);
    }

    drawLowerRightHT(x2, y2, kMinWidthT, kagekMinWidthY) {//for T design
      var poly = new Polygon();
      poly.push(x2, y2 + kagekMinWidthY);
      poly.push(x2 + kMinWidthT * 0.5, y2 - kagekMinWidthY * 4);
      poly.push(x2 + kMinWidthT * 2, y2 - kagekMinWidthY);
      poly.push(x2 + kMinWidthT * 2, y2 + kagekMinWidthY);
      this.polygons.push(poly);
    }

    drawNewGTHbox(x2m, y2, kMinWidthT, kagekMinWidthY) {//for new GTH box's left bottom corner MUKI KANKEINASHI
      var poly = new Polygon();
      poly.push(x2m, y2 - kagekMinWidthY * 5);
      poly.push(x2m - kMinWidthT * 2, y2);
      poly.push(x2m - kagekMinWidthY, y2 + kagekMinWidthY * 5);
      poly.push(x2m + kMinWidthT, y2 + kagekMinWidthY);
      poly.push(x2m, y2);
      poly.reverse();
      this.polygons.push(poly);
    }

    drawNewGTHbox_v(x2, y2, kMinWidthT, kagekMinWidthY) {
      var poly = new Polygon();
      poly.push(x2 - kMinWidthT, y2 - kagekMinWidthY * 3);
      poly.push(x2 - kMinWidthT * 2, y2);
      poly.push(x2 - kagekMinWidthY, y2 + kagekMinWidthY * 5);
      poly.push(x2 + kMinWidthT, y2 + kagekMinWidthY);
      poly.reverse();
      this.polygons.push(poly);
    }
    
    drawTailCircle_tan(tailX, tailY, dir, r, tan1, tan2) {
      //draw a (semi)circle on the tail of the line to (tailX, tailY)
      var poly = new Polygon();
      const vec1 = vector_to_len(tan1, r*bez_cir*0.74);
      const vec2 = vector_to_len(tan2, r*bez_cir*0.78);
      let p = new PointMaker(tailX, tailY, dir, r);
      let [x1, y1] = p.vec(0, -1);
      let [x2, y2] = p.vec(0, 1);
      poly.push(x1, y1);
      poly.push(x1 + vec1[0], y1 + vec1[1], 2);
      poly.push2(p.vec(0.94,-bez_cir*1.09), 2);
      poly.push2(p.vec(0.94,0));
      poly.push2(p.vec(0.94,+bez_cir*1.09), 2);
      poly.push(x2 + vec2[0], y2 + vec2[1], 2);
      poly.push(x2, y2);
      this.polygons.push(poly);
    }

    drawTailCircle(tailX, tailY, dir, r) {
      //draw a (semi)circle on the tail of the line to (tailX, tailY)
      var poly = new Polygon();
      var p = new PointMaker(tailX, tailY, dir, r);
      poly.push2(p.vec(0, -1));
      poly.push2(p.vec(bez_cir, -1), 2);
      poly.push2(p.vec(1, -bez_cir), 2);
      poly.push2(p.vec(1, 0));
      poly.push2(p.vec(1, bez_cir), 2);
      poly.push2(p.vec(bez_cir, 1), 2);
      poly.push2(p.vec(0, 1));
      this.polygons.push(poly);
    }

    drawCircle_bend_pos(x2, y2, dir, kMinWidthT2) {
      var poly = new Polygon();
      var p = new PointMaker(x2, y2, dir, kMinWidthT2);
      poly.push2(p.vec(0, -1));
      poly.push2(p.vec(1.5, -1), 2);
      poly.push2(p.vec(0.9,  1), 2);
      poly.push2(p.vec(-1,  1));
      this.polygons.push(poly);
    }

    drawCircle_bend_neg(x2, y2, dir, kMinWidthT2) {
      var poly = new Polygon();
      var p = new PointMaker(x2, y2, dir, kMinWidthT2);
      poly.push2(p.vec(0, 1));
      poly.push2(p.vec(1.5, 1), 2);
      poly.push2(p.vec(0.9,  -1), 2);
      poly.push2(p.vec(-1,  -1));
      poly.reverse();
      this.polygons.push(poly);
    }

    drawL2RSweepEnd(x2, y2, dir, kMinWidthT, kagekL2RDfatten) {
      var type = (Math.atan2(Math.abs(dir.sin), Math.abs(dir.cos)) / Math.PI * 2 - 0.6);
      if (type > 0) {
        type = type * 8;
      } else {
        type = type * 3;
      }
      const pm = (type < 0) ? -1 : 1;
      var p = new PointMaker(x2, y2, dir);
      var poly = new Polygon();
      poly.push2(p.vec(0, kMinWidthT * kagekL2RDfatten));
      poly.push2(p.vec(0, -kMinWidthT * kagekL2RDfatten));
      poly.push2(p.vec(kMinWidthT * kagekL2RDfatten * Math.abs(type), kMinWidthT * kagekL2RDfatten * pm));
      this.polygons.push(poly);
    }

    drawTurnUpwards_pos(x2, y2, kMinWidthT, length_param, dir) {
      var poly = new Polygon();
      var p = new PointMaker(x2, y2, dir);
      poly.push2(p.vec(kMinWidthT*0.7,  - kMinWidthT*0.7));
      poly.push2(p.vec(0, - kMinWidthT),2);
      poly.push2(p.vec(kMinWidthT*0.3, - kMinWidthT - length_param/2), 2);
      poly.push2(p.vec(kMinWidthT*0.5, - kMinWidthT - length_param));
      poly.push2(p.vec(0,  - kMinWidthT - length_param));
      poly.push2(p.vec( - kMinWidthT*0.6, - kMinWidthT - length_param/4), 2);
      poly.push2(p.vec( - kMinWidthT*1.8, - kMinWidthT), 2);
      poly.push2(p.vec( - kMinWidthT*2.2, - kMinWidthT));
      
      poly.reverse(); // for fill-rule
      this.polygons.push(poly);
    }

    drawTurnUpwards_neg(x2, y2, kMinWidthT, length_param, dir) {
      var poly = new Polygon();
      var p = new PointMaker(x2, y2, dir);
      poly.push2(p.vec(kMinWidthT*0.7, kMinWidthT*0.7));
      poly.push2(p.vec(0, kMinWidthT),2);
      poly.push2(p.vec(kMinWidthT*0.3, kMinWidthT + length_param/2), 2);
      poly.push2(p.vec(kMinWidthT*0.5, kMinWidthT + length_param));
      poly.push2(p.vec(0, kMinWidthT + length_param));
      poly.push2(p.vec( - kMinWidthT*0.6, kMinWidthT + length_param/4), 2);
      poly.push2(p.vec( - kMinWidthT*1.8, kMinWidthT), 2);
      poly.push2(p.vec( - kMinWidthT*2.2, kMinWidthT));
      this.polygons.push(poly);
    }
    
    drawTurnLeft(x2, y2, kMinWidthT, length_param) {
      var poly = new Polygon();
      poly.push(x2, y2 - kMinWidthT);
      poly.push(x2 - length_param, y2 - kMinWidthT);
      poly.push(x2 - length_param, y2 - kMinWidthT * 0.5);
      poly.push(x2 - kMinWidthT*0.7, y2 - kMinWidthT * 0.3,2);
      poly.push(x2, y2 + kMinWidthT * 0.3, 2);
      poly.push(x2, y2 + kMinWidthT);
      poly.reverse();
      this.polygons.push(poly);
    }

    drawUroko_h(x2, y2, kagekMinWidthY, urokoX, urokoY) {
      var poly = new Polygon();
      poly.push(x2, y2 - kagekMinWidthY);
      poly.push(x2 - urokoX, y2);
      poly.push(x2 - urokoX / 2, y2 - urokoY);
      this.polygons.push(poly);
    }

    drawUroko(x2, y2, dir, kagekMinWidthY, urokoX, urokoY) {
      var poly = new Polygon();
      poly.push(x2 + dir.sin * kagekMinWidthY, y2 - dir.cos * kagekMinWidthY);
      poly.push(x2 - dir.cos * urokoX, y2 - dir.sin * urokoX);
      poly.push(x2 - dir.cos * urokoX / 2 + dir.sin * urokoY, y2 - dir.sin * urokoX / 2 - dir.cos * urokoY);
      this.polygons.push(poly);
    }

    drawLine(x1, y1, x2, y2, halfWidth) {
      //draw a line(rectangle).
      var poly = new Polygon;
      const dir = get_dir(x2-x1, y2-y1);
      let p = new PointMaker(x1, y1, dir);
      poly.push2(p.vec(0, halfWidth));
      poly.push2(p.vec(0, -halfWidth));
      p.setpos(x2, y2);
      poly.push2(p.vec(0, -halfWidth));
      poly.push2(p.vec(0, halfWidth));
      this.polygons.push(poly);
    }

    drawOffsetLine(x1, y1, x2, y2, halfWidth, off_left1, off_right1, off_left2, off_right2) {
      //Draw a line(rectangle), and adjust the rectangle to trapezoid.
      //off_left1 means how much the start point of the left side (seeing from (x1, y1)) of the rectangle is stretched.
      //Note that the positive Y-axis goes in the downwards.
      //off_left2 deals with the end point.  The right side is managed similarly.
      //The generated polygon will be clockwise.
      var poly = new Polygon;
      const dir = get_dir(x2-x1, y2-y1);
      let p = new PointMaker(x1, y1, dir);
      poly.push2(p.vec(-off_right1, halfWidth));
      poly.push2(p.vec(-off_left1, -halfWidth));
      p.setpos(x2, y2);
      poly.push2(p.vec(off_left2, -halfWidth));
      poly.push2(p.vec(off_right2, halfWidth));
      this.polygons.push(poly);
    }

    drawCBezier(x1, y1, sx1, sy1, sx2, sy2, x2, y2, width_func, width_func_d, curve_step) {
      let [bez1, bez2] = Bezier.cBezier(x1, y1, sx1, sy1, sx2, sy2, x2, y2, width_func, width_func_d, curve_step);
      var poly = Bezier.bez_to_poly(bez1);
      poly.concat(Bezier.bez_to_poly(bez2));
      this.polygons.push(poly);
    }
    drawQBezier(x1, y1, sx, sy, x2, y2, width_func, width_func_d, curve_step) {
      let [bez1, bez2] = Bezier.qBezier(x1, y1, sx, sy, x2, y2, width_func, width_func_d, curve_step);
      var poly = Bezier.bez_to_poly(bez1);
      poly.concat(Bezier.bez_to_poly(bez2));
      this.polygons.push(poly);
    }
    drawQBezier2(x1, y1, sx, sy, x2, y2, width_func, width_func_d, curve_step) {
      let [bez1, bez2] = Bezier.qBezier2(x1, y1, sx, sy, x2, y2, width_func, width_func_d, curve_step);
      var poly = Bezier.bez_to_poly(bez1);
      poly.concat(Bezier.bez_to_poly(bez2));
      this.polygons.push(poly);
    }
  }

  class Gothic{
    constructor(size) {
      this.kRate = 50;
      if (size == 1) {
        this.kWidth = 3;
        this.kKakato = 1.8;
        this.kMage = 6;
      } else {
        this.kWidth = 5;
        this.kKakato = 3;
        this.kMage = 10;
      }
    }
    getPolygons(glyphData) {
      var cv = new FontCanvas();
      for (let glyph of glyphData) {
        this.drawStroke(cv, glyph);
      }
      return cv.getPolygons();
    }
    drawStroke(cv, s){ // gothic
      const a1 = s[0];
      const a2 = s[1];
      const a3 = s[2];
      const x1 = s[3];
      const y1 = s[4];
      const x2 = s[5];
      const y2 = s[6];
      const x3 = s[7];
      const y3 = s[8];
      const x4 = s[9];
      const y4 = s[10];
      1000 / this.kRate;
      switch(a1 % 100){
      case 0:
        break;
      case 1:
        if(a3 == 4){
          let [tx1, ty1] = get_extended_dest(x2, y2, x1, y1, -this.kMage);
          this.gothicDrawLine(x1, y1, tx1, ty1, a2, 1, cv);
          this.gothicDrawCurve(tx1, ty1, x2, y2, x2 - this.kMage * 2, y2 - this.kMage * 0.5, a1, a2, cv);
        }
        else {
          this.gothicDrawLine(x1, y1, x2, y2, a2, a3, cv);
        }
        break;
      case 2:
      case 12:{
        if(a3 == 4){
          let [tx1, ty1] = get_extended_dest_wrong(x3, y3, x2, y2, -this.kMage);
          this.gothicDrawCurve(x1, y1, x2, y2, tx1, ty1, a1, a2, cv);
          this.gothicDrawCurve(tx1, ty1, x3, y3, x3 - this.kMage * 2, y3 - this.kMage * 0.5, a1, a2, cv);
        }
        else if(a3 == 5){
          const tx1 = x3 + this.kMage;
          const ty1 = y3;
          const tx2 = tx1 + this.kMage * 0.5;
          const ty2 = y3 - this.kMage * 2;
          this.gothicDrawCurve(x1, y1, x2, y2, x3, y3, a1, a2, cv);
          this.gothicDrawCurve(x3, y3, tx1, ty1, tx2, ty2, a1, a2, cv);
        }
        else {
          this.gothicDrawCurve(x1, y1, x2, y2, x3, y3, a1, a2, cv);
        }
        break;
      }
      case 3:{
        let [tx1, ty1] = get_extended_dest(x2, y2, x1, y1, -this.kMage);
          let [tx2, ty2] = get_extended_dest(x2, y2, x3, y3, -this.kMage);
          
        if(a3 == 5){
          const tx3 = x3 - this.kMage;
          const ty3 = y3;
          const tx4 = x3 + this.kMage * 0.5;
          const ty4 = y3 - this.kMage * 2;
          this.gothicDrawLine(x1, y1, tx1, ty1, a2, 1, cv);
          this.gothicDrawCurve(tx1, ty1, x2, y2, tx2, ty2, a1, a2, cv);
          this.gothicDrawLine(tx2, ty2, tx3, ty3, 1, 1, cv);
          this.gothicDrawCurve(tx3, ty3, x3, y3, tx4, ty4, a1, a2, cv);
        }
        else {
          this.gothicDrawLine(x1, y1, tx1, ty1, a2, 1, cv);
          this.gothicDrawCurve(tx1, ty1, x2, y2, tx2, ty2, a1, a2, cv);
          this.gothicDrawLine(tx2, ty2, x3, y3, 1, a3, cv);
        }
        break;
      }
      case 6:
        if(a3 == 5){
          const tx1 = x4 - this.kMage;
          const ty1 = y4;
          const tx2 = x4 + this.kMage * 0.5;
          const ty2 = y4 - this.kMage * 2;
          cv.drawCBezier(x1, y1, x2, y2, x3, y3, tx1, ty1, (t) => { return this.kWidth; }, t => 0, 1000 / this.kRate);
          this.gothicDrawCurve(tx1, ty1, x4, y4, tx2, ty2, a1, a2, cv);
        }
        else {
          cv.drawCBezier(x1, y1, x2, y2, x3, y3, x4, y4,  (t) => { return this.kWidth; }, t => 0, 1000 / this.kRate);
        }
        break;
      case 7:
        this.gothicDrawLine(x1, y1, x2, y2, a2, 1, cv);
        this.gothicDrawCurve(x2, y2, x3, y3, x4, y4, a1, a2, cv);
        break;
      }
    }

    gothicDrawCurve(x1, y1, x2, y2, x3, y3, ta1, ta2, cv) {
      cv.drawQBezier(x1, y1, x2, y2, x3, y3, (t) => { return this.kWidth; }, t => 0, 1000 / this.kRate);
    }

    gothicDrawLine(tx1, ty1, tx2, ty2, ta1, ta2, cv) {
      var x1 = tx1;
      var y1 = ty1;
      var x2 = tx2;
      var y2 = ty2;
      if (ta1 % 10 == 2) {
        let [x1ext, y1ext] = get_extended_dest(tx1, ty1, tx2, ty2, this.kWidth);
        x1 = x1ext; y1 = y1ext;
      } else if (ta1 % 10 == 3) {
        let [x1ext, y1ext] = get_extended_dest(tx1, ty1, tx2, ty2, this.kWidth * this.kKakato);
        x1 = x1ext; y1 = y1ext;
      }
      if (ta2 % 10 == 2) {
        let [x2ext, y2ext] = get_extended_dest(tx2, ty2, tx1, ty1, this.kWidth);
        x2 = x2ext; y2 = y2ext;
      } else if (ta2 % 10 == 3) {
        let [x2ext, y2ext] = get_extended_dest(tx2, ty2, tx1, ty1, this.kWidth * this.kKakato);
        x2 = x2ext; y2 = y2ext;
      }
      cv.drawLine(x1, y1, x2, y2, this.kWidth);
    }
  }

  const STROKETYPE = {
      STRAIGHT : 1,//直線
      CURVE : 2,//曲線 quadratic bezier curve
      BENDING : 3,//折れ used in "札" "己"
      BENDING_ROUND : 4,//used for "乙"
      BEZIER : 6,//複曲線 cubic bezier curve
      VCURVE : 7,//vertical line and curve. used in the leftmost stroke of "月".
                 //although the angle of line can be chosen arbitrarily, only vertical line is expected.
      REFERENCE : 99,
  };
  const STARTTYPE = {
      OPEN : 0,//simple lines like "三" or "川" (two strokes on the right side)
               //also used in the left stroke of "人"
      CONNECTING_H : 2,//horizontal strokes connecting to other strokes.  used in the center strokes of "日".
      UPPER_LEFT_CORNER : 12,//the starting point is at the upper left corner.  usually used for vertical lines, like the leftmost stroke of "日".
      UPPER_RIGHT_CORNER : 22,//the starting point is at the upper right corner.  usually used for vertical lines, like the rightmost stroke of "日".
      CONNECTING_V : 32,//vertical strokes connecting to other strokes.  used in the center strokes of "工".
      THIN : 7//used in the right stroke of "人"
  };
  const ENDTYPE = {
      OPEN : 0,//simple lines like "三" or "川" (two strokes on the right side)
               //also used in the right stroke of "人"(L2R sweep)
      CONNECTING_H : 2,
      CONNECTING_V : 32,///vertical strokes connecting to other strokes.  used in the center strokes of "工".
      LOWER_LEFT_CORNER : 13,
      LOWER_RIGHT_CORNER : 23,
      LOWER_LEFT_ZH_OLD : 313,//for characters used in China.
      LOWER_LEFT_ZH_NEW : 413,//for characters used in China.
      
      TURN_LEFT : 4,//adds a short line to the left.  used in the rightmost stroke of "月".
      LOWER_RIGHT_HT : 24,//for characters used in China.
      TURN_UPWARDS : 5,//adds a short upward line.  used in the rightmost stroke of "札" or "風".
      LEFT_SWEEP : 7, //used in the left stroke of "人".
      STOP : 8,//used in the rightmost stroke of "小" or lower four dots of "魚".
  };

  // Reference : http://www.cam.hi-ho.ne.jp/strong_warriors/teacher/chapter0{4,5}.html

  function point(x, y){
    this.x = x;
    this.y = y;
  }

  function getCrossPoint(x11, y11, x12, y12, x21, y21, x22, y22){ // point
    var a1 = y12 - y11;
    var b1 = x11 - x12;
    var c1 = -1 * a1 * x11 - b1 * y11;
    var a2 = y22 - y21;
    var b2 = x21 - x22;
    var c2 = -1 * a2 * x21 - b2 * y21;
    
    var temp = b1 * a2 - b2 * a1;
    if(temp == 0){ // parallel
      return false;
    }
    return new point((c1 * b2 - c2 * b1) / temp, (a1 * c2 - a2 * c1) / temp);
  }

  function isCross(x11, y11, x12, y12, x21, y21, x22, y22){ // boolean
    var temp = getCrossPoint(x11, y11, x12, y12, x21, y21, x22, y22);
    if(!temp){ return false; }
    if(x11 < x12 && (temp.x < x11 || x12 < temp.x) ||
       x11 > x12 && (temp.x < x12 || x11 < temp.x) ||
       y11 < y12 && (temp.y < y11 || y12 < temp.y) ||
       y11 > y12 && (temp.y < y12 || y11 < temp.y)
       ){
      return false;
    }
    if(x21 < x22 && (temp.x < x21 || x22 < temp.x) ||
       x21 > x22 && (temp.x < x22 || x21 < temp.x) ||
       y21 < y22 && (temp.y < y21 || y22 < temp.y) ||
       y21 > y22 && (temp.y < y22 || y21 < temp.y)
       ){
      return false;
    }
    return true;
  }

  function isCrossBox(x1, y1, x2, y2, bx1, by1, bx2, by2){ // boolean
    if(isCross(x1, y1, x2, y2, bx1, by1, bx2, by1)){ return true; }
    if(isCross(x1, y1, x2, y2, bx2, by1, bx2, by2)){ return true; }
    if(isCross(x1, y1, x2, y2, bx1, by2, bx2, by2)){ return true; }
    if(isCross(x1, y1, x2, y2, bx1, by1, bx1, by2)){ return true; }
    return false;
  }

  function isCrossBoxWithOthers(strokesArray, i, bx1, by1, bx2, by2){ // boolean
    for(var j = 0; j < strokesArray.length; j++){
      if(i == j){ continue; }
      switch(strokesArray[j][0]){
      case 0:
      case 8:
      case 9:
        break;
      case 6:
      case 7:
        if(isCrossBox(strokesArray[j][7],
                      strokesArray[j][8],
                      strokesArray[j][9],
                      strokesArray[j][10],
                      bx1, by1, bx2, by2)){
          return true;
        }
      case 2:
      case 12:
      case 3:
        if(isCrossBox(strokesArray[j][5],
                      strokesArray[j][6],
                      strokesArray[j][7],
                      strokesArray[j][8],
                      bx1, by1, bx2, by2)){
          return true;
        }
      default:
        if(isCrossBox(strokesArray[j][3],
                      strokesArray[j][4],
                      strokesArray[j][5],
                      strokesArray[j][6],
                      bx1, by1, bx2, by2)){
          return true;
        }
      }
    }
    return false;
  }

  function isCrossWithOthers(strokesArray, i, bx1, by1, bx2, by2){ // boolean
    for(var j = 0; j < strokesArray.length; j++){
      if(i == j){ continue; }
      switch(strokesArray[j][0]){
      case 0:
      case 8:
      case 9:
        break;
      case 6:
      case 7:
        if(isCross(strokesArray[j][7],
                   strokesArray[j][8],
                   strokesArray[j][9],
                   strokesArray[j][10],
                   bx1, by1, bx2, by2)){
          return true;
        }
      case 2:
      case 12:
      case 3:
        if(isCross(strokesArray[j][5],
                   strokesArray[j][6],
                   strokesArray[j][7],
                   strokesArray[j][8],
                   bx1, by1, bx2, by2)){
          return true;
        }
      default:
        if(isCross(strokesArray[j][3],
                   strokesArray[j][4],
                   strokesArray[j][5],
                   strokesArray[j][6],
                   bx1, by1, bx2, by2)){
          return true;
        }
      }
    }
    return false;
  }

  class Mincho {
    constructor(size) {
      this.kRate = 50;
      if (size == 1) {
        this.kMinWidthY = 1.2;
        this.kMinWidthT = 3.6;
        this.kWidth = 3;
        this.kKakato = 1.8;
        this.kL2RDfatten = 1.1;
        this.kMage = 6;

        this.kAdjustKakatoL = ([8, 5, 3, 1]); // for KAKATO adjustment 000,100,200,300,400
        this.kAdjustKakatoR = ([4, 3, 2, 1]); // for KAKATO adjustment 000,100,200,300
        this.kAdjustKakatoRangeX = 12; // check area width
        this.kAdjustKakatoRangeY = ([1, 11, 14, 18]); // 3 steps of checking
        this.kAdjustKakatoStep = 3; // number of steps

        this.kAdjustUrokoX = ([14, 12, 9, 7]); // for UROKO adjustment 000,100,200,300
        this.kAdjustUrokoY = ([7, 6, 5, 4]); // for UROKO adjustment 000,100,200,300
        this.kAdjustUrokoLength = ([13, 21, 30]); // length for checking
        this.kAdjustUrokoLengthStep = 3; // number of steps
        this.kAdjustUrokoLine = ([13, 15, 18]); // check for crossing. corresponds to length

        this.kAdjustUroko2Step = 3;
        this.kAdjustUroko2Length = 25;

        this.kAdjustTateStep = 4;
        this.kAdjustMageStep = 5;
      } else if (size == 3) {
        this.kMinWidthY = 3;
        this.kMinWidthT = 8;
        this.kWidth = 6;
        this.kKakato = 4;
        this.kL2RDfatten = 1.1;
        this.kMage = 14;

        this.kAdjustKakatoL = ([20, 13, 7, 3]); // for KAKATO adjustment 000,100,200,300,400
        this.kAdjustKakatoR = ([12, 9, 6, 3]); // for KAKATO adjustment 000,100,200,300
        this.kAdjustKakatoRangeX = 26; // check area width
        this.kAdjustKakatoRangeY = ([2, 26, 30, 40]); // 3 steps of checking
        this.kAdjustKakatoStep = 3; // number of steps

        this.kAdjustUrokoX = ([30, 25, 20, 15]); // for UROKO adjustment 000,100,200,300
        this.kAdjustUrokoY = ([15, 14, 13, 12]); // for UROKO adjustment 000,100,200,300
        this.kAdjustUrokoLength = ([29, 40, 62]); // length for checking
        this.kAdjustUrokoLengthStep = 3; // number of steps
        this.kAdjustUrokoLine = ([29, 34, 39]); // check for crossing. corresponds to length

        this.kAdjustUroko2Step = 3;
        this.kAdjustUroko2Length = 50;

        this.kAdjustTateStep = 4;
        this.kAdjustMageStep = 5;
      } else if (size > 1) {
        this.kMinWidthY = size;
        this.kMinWidthT = size * 2.6;
        this.kWidth = size * 2.2;
        this.kKakato = size * 1.2 + 0.6;
        this.kL2RDfatten = 1.2;
        this.kMage = size * 4 + 2;

        this.kAdjustKakatoL = ([size * 4 + 1, size * 3+0.75, size * 2 + 0.5, size * 1 + 0.25]); // for KAKATO adjustment 000,100,200,300,400
        this.kAdjustKakatoR = ([size * 3.2, size * 2.4, size * 1.6, size * 0.8]); // for KAKATO adjustment 000,100,200,300
        this.kAdjustKakatoRangeX = size * 4 + 8; // check area width
        this.kAdjustKakatoRangeY = ([size, size * 5 + 3, size * 9 + 5, size * 12 + 6]); // 3 steps of checking
        this.kAdjustKakatoStep = 3; // number of steps

        this.kAdjustUrokoX = ([size * 8 + 4, size * 7 + 3.5, size * 6 + 3, size * 5 + 2.5]); // for UROKO adjustment 000,100,200,300
        this.kAdjustUrokoY = ([size * 5 + 2, size * 4.8 + 1.5, size * 4.6 + 1, size * 4.4 + 0.5]); // for UROKO adjustment 000,100,200,300
        this.kAdjustUrokoLength = ([size * 7 + 7, size * 11 + 11, size * 15 + 15]); // length for checking
        this.kAdjustUrokoLengthStep = 3; // number of steps
        this.kAdjustUrokoLine = ([size * 7 + 7, size * 9 + 8, size * 11 + 9]); // check for crossing. corresponds to length

        this.kAdjustUroko2Step = 3;
        this.kAdjustUroko2Length = size * 20;

        this.kAdjustTateStep = 4;
        this.kAdjustMageStep = 5;
      } else {
        this.kMinWidthY = 2;
        this.kMinWidthT = 6;
        this.kWidth = 5;
        this.kKakato = 3;
        this.kL2RDfatten = 1.1;
        this.kMage = 10;

        this.kAdjustKakatoL = ([14, 9, 5, 2]); // for KAKATO adjustment 000,100,200,300,400
        this.kAdjustKakatoR = ([8, 6, 4, 2]); // for KAKATO adjustment 000,100,200,300
        this.kAdjustKakatoRangeX = 20; // check area width
        this.kAdjustKakatoRangeY = ([1, 19, 24, 30]); // 3 steps of checking
        this.kAdjustKakatoStep = 3; // number of steps

        this.kAdjustUrokoX = ([24, 20, 16, 12]); // for UROKO adjustment 000,100,200,300
        this.kAdjustUrokoY = ([12, 11, 9, 8]); // for UROKO adjustment 000,100,200,300
        this.kAdjustUrokoLength = ([22, 36, 50]); // length for checking
        this.kAdjustUrokoLengthStep = 3; // number of steps
        this.kAdjustUrokoLine = ([22, 26, 30]); // check for crossing. corresponds to length

        this.kAdjustUroko2Step = 3;
        this.kAdjustUroko2Length = 40;

        this.kAdjustTateStep = 4;
        this.kAdjustMageStep = 5;
      }
    }

    getPolygons(glyphData) {
      var cv = new FontCanvas();
      for (var i = 0; i < glyphData.length; i++) {
        var tempdata = glyphData.slice();
        tempdata.splice(i, 1);
        this.drawAdjustedStroke(cv, glyphData[i], tempdata);
      }
      return cv.getPolygons();
    }

    drawAdjustedStroke(cv, s, others) {//draw stroke on the canvas
      const a1 = s[0];
      const a2 = s[1];
      const a3 = s[2];
      const x1 = s[3];
      const y1 = s[4];
      const x2 = s[5];
      const y2 = s[6];
      const x3 = s[7];
      const y3 = s[8];
      const x4 = s[9];
      const y4 = s[10];
      if(a2>100){
        console.log("error: start type"+a2);
      }
      if(a3>100){
        console.log("error: end type"+a3);
      }
      
      const dir12 = get_dir(x2-x1, y2-y1);
      const dir23 = get_dir(x3-x2, y3-y2);
      const dir34 = get_dir(x4-x3, y4-y3);
      
      switch (a1) {
        case 0: {//rotate and flip
          if (a2 == 98) {
            cv.flip_left_right(x1, y1, x2, y2);
          } else if (a2 == 97) {
            cv.flip_up_down(x1, y1, x2, y2);
          } else if (a2 == 99 && a3 == 1) {
            cv.rotate90(x1, y1, x2, y2);
          } else if (a2 == 99 && a3 == 2) {
            cv.rotate180(x1, y1, x2, y2);
          } else if (a2 == 99 && a3 == 3) {
            cv.rotate270(x1, y1, x2, y2);
          }
          break;
        }
        case STROKETYPE.STRAIGHT: {
          const dir = get_dir(x2-x1, y2-y1);
          if (a3 == ENDTYPE.CONNECTING_H) {//usually horizontal
            cv.drawLine(x1, y1, x2, y2, this.kMinWidthY);
          } else if (a3 == ENDTYPE.OPEN && Math.abs(y2 - y1) < x2 - x1) { //horizontal or gentle slope
            const param_uroko = this.adjustUrokoParam(s, others);
            const param_uroko2 = this.adjustUroko2Param(s, others);
            cv.drawLine(x1, y1, x2, y2, this.kMinWidthY);
            if (y1 == y2) {//horizontal
              const uroko_max = Math.max(param_uroko, param_uroko2);
              cv.drawUroko_h(x2, y2, this.kMinWidthY, this.kAdjustUrokoX[uroko_max], this.kAdjustUrokoY[uroko_max]);
            } else {
              cv.drawUroko(x2, y2, dir, this.kMinWidthY, this.kAdjustUrokoX[param_uroko], this.kAdjustUrokoY[param_uroko]);
            }
          } else {//vertical or steep slope
            let poly_end = new Polygon(2);
            const param_tate = this.adjustTateParam(s, others);
            const kMinWidthT_m = this.kMinWidthT - param_tate / 2;
            //head
            let poly_start = this.getStartOfVLine(x1, y1, x2, y2, a2, kMinWidthT_m, cv);
            //tail
            switch (a3) {
              case ENDTYPE.OPEN: {
                const right2 = kMinWidthT_m / 2;
                const left2 = -kMinWidthT_m / 2;
                poly_end = this.getEndOfOffsetLine(x1, y1, x2, y2, kMinWidthT_m, right2, left2);
                break;
              }
              case ENDTYPE.TURN_LEFT: {
                let [tx1, ty1] = moved_point(x2, y2, dir12, -this.kMage);
                const width_func = (t) => { return kMinWidthT_m; };
                const new_x2 = x2 - this.kMage * (((this.kAdjustTateStep + 4) - param_tate) / (this.kAdjustTateStep + 4));
                cv.drawQBezier(tx1, ty1, x2, y2,
                  new_x2, y2, width_func, t => 0);
                const param_hane = this.adjustHaneParam(x2, y2, others);
                cv.drawTurnLeft(new_x2, y2, kMinWidthT_m, this.kWidth * 4 * Math.min(1 - param_hane / 10, Math.pow(kMinWidthT_m / this.kMinWidthT, 3)));
                poly_end = this.getEndOfLine(x1, y1, tx1, ty1, kMinWidthT_m);
                break;
              }
              case ENDTYPE.LOWER_LEFT_CORNER: {
                const param_kakato = this.adjustKakatoParam(s, others);
                const right2 = this.kAdjustKakatoL[param_kakato] + kMinWidthT_m;
                const left2 = this.kAdjustKakatoL[param_kakato];
                poly_end = this.getEndOfOffsetLine(x1, y1, x2, y2, kMinWidthT_m, right2, left2);
                break;
              }
              case ENDTYPE.LOWER_RIGHT_CORNER: {
                const param_kakato = this.adjustKakatoParam(s, others);
                const right2 = this.kAdjustKakatoR[param_kakato] + kMinWidthT_m;
                const left2 = this.kAdjustKakatoR[param_kakato];
                poly_end = this.getEndOfOffsetLine(x1, y1, x2, y2, kMinWidthT_m, right2, left2);
                break;
              }
              case ENDTYPE.LOWER_LEFT_ZH_NEW: {
                if (x1 == x2) {//vertical
                  cv.drawNewGTHbox_v(x2, y2, kMinWidthT_m, this.kMinWidthY);
                } else {
                  var m = 0;
                  if (x1 > x2 && y1 != y2) {
                    m = Math.floor((x1 - x2) / (y2 - y1) * 3);
                  }
                  cv.drawNewGTHbox(x2 + m, y2, kMinWidthT_m, this.kMinWidthY);
                }
                const right2 = kMinWidthT_m;
                const left2 = 0;
                poly_end = this.getEndOfOffsetLine(x1, y1, x2, y2, kMinWidthT_m, right2, left2);
                break;
              }
              case ENDTYPE.LOWER_LEFT_ZH_OLD: {
                const right2 = this.kAdjustKakatoL[3] + kMinWidthT_m;
                const left2 = this.kAdjustKakatoL[3];
                poly_end = this.getEndOfOffsetLine(x1, y1, x2, y2, kMinWidthT_m, right2, left2);
                break;
              }
              case ENDTYPE.CONNECTING_V: {
                if (y1 == y2) {//horizontal (error)
                  console.log("error: connecting_v at the end of the horizontal line");
                  cv.drawLine(x1, y1, x2, y2, kMinWidthT_m);
                } else if (x1 == x2) {//vertical
                  poly_end.set(0, x2 + kMinWidthT_m, y2 + this.kMinWidthY - 0.001);
                  poly_end.set(1, x2 - kMinWidthT_m, y2 + this.kMinWidthY - 0.001);
                } else {
                  const rad = Math.atan((y2 - y1) / (x2 - x1));
                  const v = (x1 > x2) ? -1 : 1;
                  poly_end.set(0, x2 + (kMinWidthT_m * v) / Math.sin(rad), y2);
                  poly_end.set(1, x2 - (kMinWidthT_m * v) / Math.sin(rad), y2);
                }
                break;
              }
              case ENDTYPE.LOWER_RIGHT_HT: {
                if (x1 == x2) {//vertical
                  cv.drawLowerRightHT_v(x2, y2, kMinWidthT_m, this.kMinWidthY);
                } else {
                  cv.drawLowerRightHT(x2, y2, kMinWidthT_m, this.kMinWidthY);
                }
                if (y1 == y2) {//horizontal (error)
                  console.log("error: connecting_v at the end of the horizontal line");
                  cv.drawLine(x1, y1, x2, y2, kMinWidthT_m);
                } else if (x1 == x2) {//vertical
                  poly_end.set(0, x2 + kMinWidthT_m, y2 + this.kMinWidthY);
                  poly_end.set(1, x2 - kMinWidthT_m, y2 + this.kMinWidthY);
                } else {
                  const rad = Math.atan((y2 - y1) / (x2 - x1));
                  const v = (x1 > x2) ? -1 : 1;
                  poly_end.set(0, x2 + (kMinWidthT_m * v) / Math.sin(rad), y2);
                  poly_end.set(1, x2 - (kMinWidthT_m * v) / Math.sin(rad), y2);
                }
                break;
              }
              default:
                throw ("error: unknown end type at the straight line");
            }
            //body
            poly_start.concat(poly_end);
            cv.addPolygon(poly_start);
          }
          break;
        }



        case STROKETYPE.CURVE: {
          //head
          if (a2 == STARTTYPE.OPEN) {
            let [x1ext, y1ext] = moved_point(x1, y1, dir12, 1 * this.kMinWidthY * 0.5);
            if (y1ext <= y3) { //from up to bottom
              cv.drawOpenBegin_curve_down(x1ext, y1ext, dir12, this.kMinWidthT, this.kMinWidthY);
            }
            else { //from bottom to up
              cv.drawOpenBegin_curve_up(x1ext, y1ext, dir12, this.kMinWidthT, this.kMinWidthY);
            }
          } else if (a2 == STARTTYPE.UPPER_RIGHT_CORNER) {
            cv.drawUpperRightCorner(x1, y1, this.kMinWidthT, this.kMinWidthY, this.kWidth);
          } else if (a2 == STARTTYPE.UPPER_LEFT_CORNER) {
            let [x1ext, y1ext] = moved_point(x1, y1, dir12, -this.kMinWidthY);
            cv.drawUpperLeftCorner(x1ext, y1ext, dir12, this.kMinWidthT);
          }
          //body
          const a2temp = (a2 == STARTTYPE.CONNECTING_V && this.adjustKirikuchiParam(s, others)) ? 100 + a2 : a2;
          let [tan1, tan2] = this.minchoDrawCurve(x1, y1, x2, y2, x3, y3, a2temp, a3, cv);
          //tail
          switch (a3) {
            case ENDTYPE.TURN_LEFT: {
              let [tx1, ty1] = moved_point(x3, y3, dir23, -this.kMage);
              const param_hane = this.adjustHaneParam(x3, y3, others);
              const width_func = (t) => { return this.kMinWidthT; };
              cv.drawQBezier(tx1, ty1, x3, y3, x3 - this.kMage, y3, width_func, t => 0);
              cv.drawTurnLeft(x3 - this.kMage, y3, this.kMinWidthT, this.kWidth * 4 * Math.min(1 - param_hane / 10, 1));
              break;
            }
            case ENDTYPE.TURN_UPWARDS: {
              cv.drawTailCircle(x3, y3, dir23, this.kMinWidthT);
              cv.drawTurnUpwards_pos(x3, y3, this.kMinWidthT, this.kWidth*5, (y1<y3)?DIR_POSX:DIR_NEGX);
              break;
            }
            case ENDTYPE.STOP: {
              let [x3ex, y3ex] = moved_point(x3, y3, dir23, -1 * this.kMinWidthT * 0.52);
              cv.drawTailCircle_tan(x3ex, y3ex, dir23, this.kMinWidthT*1.1, tan1, tan2);
              break;
            }
            default: {
              if (a2 == STARTTYPE.THIN && a3 == ENDTYPE.OPEN) {
                cv.drawL2RSweepEnd(x3, y3, dir23, this.kMinWidthT, this.kL2RDfatten);
              }
              break;
            }
          }
          break;
        }



        case STROKETYPE.BENDING: {
          //first line
          const param_tate = this.adjustTateParam(s, others);
          const param_mage = this.adjustMageParam(s, others) / 2;
          const kMinWidthT_m = this.kMinWidthT - param_tate / 2;
          const kMinWidthT_mage = (this.kMinWidthT - param_mage / 2)*0.9;
          let [tx1, ty1] = moved_point(x2, y2, dir12, -this.kMage);
          let [tx2, ty2] = moved_point(x2, y2, dir23, this.kMage);
          {
            let poly_start = this.getStartOfVLine(x1, y1, x2, y2, a2, kMinWidthT_m, cv);
            let poly_end = this.getEndOfLine(x1, y1, tx1, ty1, kMinWidthT_m);
            poly_start.concat(poly_end);
            cv.addPolygon(poly_start);
          }
          //curve
          const width_func = function (t) {
            return kMinWidthT_mage * t + kMinWidthT_m * (1 - t);
          };
          cv.drawQBezier(tx1, ty1, x2, y2, tx2, ty2, width_func, t => 0);



          //last line
          if (tx2 < x3) {
            cv.drawOffsetLine(tx2, ty2, x3, y3, kMinWidthT_mage, 0, 0, 0, -kMinWidthT_mage);
          } else {
            cv.drawOffsetLine(tx2, ty2, x3, y3, kMinWidthT_mage, 0, 0, -kMinWidthT_mage, 0);
          }
          if (y2 == y3) {
            if (tx2 < x3) {
              cv.drawCircle_bend_pos(x3, y3, DIR_POSX, kMinWidthT_mage);
            } else {
              cv.drawCircle_bend_neg(x3, y3, DIR_NEGX, kMinWidthT_mage);
            }
            if (a3 == ENDTYPE.TURN_UPWARDS) {
              if (tx2 < x3) {
                cv.drawTurnUpwards_pos(x3, y3, kMinWidthT_mage, this.kWidth * (4 * (1 - param_mage / this.kAdjustMageStep) + 1), DIR_POSX);
              } else {
                cv.drawTurnUpwards_neg(x3, y3, kMinWidthT_mage, this.kWidth * (4 * (1 - param_mage / this.kAdjustMageStep) + 1), DIR_NEGX);
              }
            }
          } else {
            const dir = get_dir(x3-x2, y3-y2);
            if (tx2 < x3) {
              cv.drawCircle_bend_pos(x3, y3, dir, kMinWidthT_mage);
            } else {
              cv.drawCircle_bend_neg(x3, y3, dir, kMinWidthT_mage);
            }
            if (a3 == ENDTYPE.TURN_UPWARDS) {
              if (tx2 < x3) {
                cv.drawTurnUpwards_pos(x3, y3, kMinWidthT_mage, this.kWidth*5, dir);
              } else {
                cv.drawTurnUpwards_neg(x3, y3, kMinWidthT_mage, this.kWidth*5, dir);
              }
            }
          }
          break;
        }
        case 12: {
          throw "error: unknown stroketype 12";
        }
        case STROKETYPE.BENDING_ROUND: {
          var rate = 6;
          if ((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) < 14400) { // smaller than 120 x 120
            rate = Math.sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2)) / 120 * 6;
          }
          let [tx1, ty1] = moved_point(x2, y2, dir12, -this.kMage * rate);
          let [tx2, ty2] = moved_point(x2, y2, dir23,  this.kMage * rate);
          //first line
          let poly_start = this.getStartOfVLine(x1, y1, x2, y2, a2, this.kMinWidthT, cv);
          let poly_end = this.getEndOfLine(x1, y1, tx1, ty1, this.kMinWidthT);
          poly_start.concat(poly_end);
          cv.addPolygon(poly_start);
          //curve
          const width_func = (t) => { return this.kMinWidthT; };
          cv.drawQBezier(tx1, ty1, x2, y2, tx2, ty2, width_func, t => 0);


          //last line
          if (tx2 < x3) {
            cv.drawOffsetLine(tx2, ty2, x3, y3, this.kMinWidthT, 0, 0, 0, -this.kMinWidthT);
          } else {
            cv.drawOffsetLine(tx2, ty2, x3, y3, this.kMinWidthT, 0, 0, -this.kMinWidthT, 0);
          }
          if (y2 == y3) {
            if (tx2 < x3) {
              cv.drawCircle_bend_pos(x3, y3, DIR_POSX, this.kMinWidthT);
            } else {
              cv.drawCircle_bend_neg(x3, y3, DIR_NEGX, this.kMinWidthT);
            }
            if (a3 == ENDTYPE.TURN_UPWARDS) {
              if (tx2 < x3) {
                cv.drawTurnUpwards_pos(x3, y3, this.kMinWidthT, this.kWidth * 5, DIR_POSX);
              } else {
                cv.drawTurnUpwards_neg(x3, y3, this.kMinWidthT, this.kWidth * 5, DIR_NEGX);
              }
            }
          } else {
            const dir = get_dir(x3-x2, y3-y2);
            if (tx2 < x3) {
              cv.drawCircle_bend_pos(x3, y3, dir, this.kMinWidthT);
            } else {
              cv.drawCircle_bend_neg(x3, y3, dir, this.kMinWidthT);
            }
            if (a3 == ENDTYPE.TURN_UPWARDS) {
              if (tx2 < x3) {
                cv.drawTurnUpwards_pos(x3, y3, this.kMinWidthT, this.kWidth*5, dir);
              } else {
                cv.drawTurnUpwards_neg(x3, y3, this.kMinWidthT, this.kWidth*5, dir);
              }
            }
          }
          break;
        }



        case STROKETYPE.BEZIER: {
          //head
          if (a2 == STARTTYPE.OPEN) {
            let [x1ext, y1ext] = moved_point(x1, y1, dir12, 1 * this.kMinWidthY * 0.5);

            if (y1ext <= y4) { //from up to bottom
              cv.drawOpenBegin_curve_down(x1ext, y1ext, dir12, this.kMinWidthT, this.kMinWidthY);
            }
            else { //from bottom to up
              cv.drawOpenBegin_curve_up(x1ext, y1ext, dir12, this.kMinWidthT, this.kMinWidthY);
            }
          } else if (a2 == STARTTYPE.UPPER_RIGHT_CORNER) {
            cv.drawUpperRightCorner(x1, y1, this.kMinWidthT, this.kMinWidthY, this.kWidth);
          } else if (a2 == STARTTYPE.UPPER_LEFT_CORNER) {
            let [x1ext, y1ext] = moved_point(x1, y1, dir12, -this.kMinWidthY);
            cv.drawUpperLeftCorner(x1ext, y1ext, dir12, this.kMinWidthT);
          }
          //body
          let [tan1, tan2] = this.minchoDrawBezier(x1, y1, x2, y2, x3, y3, x4, y4, a2, a3, cv);
          //tail
          switch (a3) {
            case ENDTYPE.TURN_LEFT:
              let [tx1, ty1] = moved_point(x4, y4, dir34, -this.kMage);
              const width_func = (t) => { return this.kMinWidthT; };
              cv.drawQBezier(tx1, ty1, x4, y4, x4 - this.kMage, y4, width_func, t => 0);
              const param_hane = this.adjustHaneParam(x4, y4, others);
              cv.drawTurnLeft(x4 - this.kMage, y4, this.kMinWidthT, this.kWidth * 4 * Math.min(1 - param_hane / 10, 1));
              break;
            case ENDTYPE.TURN_UPWARDS:
              cv.drawTailCircle(x4, y4, dir34, this.kMinWidthT);
              cv.drawTurnUpwards_pos(x4, y4, this.kMinWidthT, this.kWidth*5, (y1<y4)?DIR_POSX:DIR_NEGX);
              break;
            case ENDTYPE.STOP:
              let [x4ex, y4ex] = moved_point(x4, y4, dir34, -this.kMinWidthT * 0.52);
              cv.drawTailCircle_tan(x4ex, y4ex, dir34, this.kMinWidthT*1.1, tan1, tan2);
              break;
            default:
              if (a2 == STARTTYPE.THIN && a3 == ENDTYPE.OPEN) {
                cv.drawL2RSweepEnd(x4, y4, dir34, this.kMinWidthT, this.kL2RDfatten);
              }
              break;
          }
          break;
        }



        case STROKETYPE.VCURVE: {
          const param_tate = this.adjustTateParam(s, others);
          const kMinWidthT_m = this.kMinWidthT - param_tate / 2;
          //straight
          let poly_start = this.getStartOfVLine(x1, y1, x2, y2, a2, kMinWidthT_m, cv);
          let poly_end = this.getEndOfLine(x1, y1, x2, y2, kMinWidthT_m);
          poly_start.concat(poly_end);
          cv.addPolygon(poly_start);
          //curve
          const width_func = function (t) {
            //const deltad = Math.pow(1.0-t,0.7)*0.8+0.2;
            const deltad = (1 - Math.pow(t, 1.8)) * 0.85 + 0.15;
            return deltad * kMinWidthT_m;
          };
          cv.drawQBezier(x2, y2, x3, y3, x4, y4, width_func, t => -1.8 * Math.pow(t, 0.8) * 0.85 * kMinWidthT_m);
          break;
        }
        case 9: // may not be exist ... no need
          //kageCanvas[y1][x1] = 0;
          //kageCanvas[y2][x2] = 0;
          break;
        default:
          throw "error: unknown stroke "+s;
      }
    }

    minchoDrawCurve(x1pre, y1pre, sx, sy, x2pre, y2pre, a1, a2, cv) {
      var delta;
      switch (a1) {
        case STARTTYPE.OPEN:
        case STARTTYPE.THIN:
          delta = -1 * this.kMinWidthY * 0.5;
          break;
        case STARTTYPE.UPPER_LEFT_CORNER:
          //case 32:
          delta = this.kMinWidthY;
          break;
        default:
          delta = 0;
          break;
      }
      let [x1, y1] = get_extended_dest(x1pre, y1pre, sx, sy, delta);

      switch (a2) {
        case ENDTYPE.STOP: // get shorten for tail's circle
          delta = -1 * this.kMinWidthT * 0.52;
          break;
        case ENDTYPE.TURN_LEFT:
          delta = -this.kMage;
          break;
        default:
          delta = 0;
          break;
      }
      let [x2, y2] = get_extended_dest(x2pre, y2pre, sx, sy, delta);

      var width_func;
      var width_func_d;
      let bez1, bez2;
      
      if (a1 == STARTTYPE.THIN && a2 == ENDTYPE.STOP) { //stop
        //const slant_cos = 
        width_func = t => widfun_stop(t, x1, y1, x2, y2, this.kMinWidthT);
        width_func_d = t => widfun_stop_d(t, x1, y1, x2, y2, this.kMinWidthT);

        [bez1, bez2] = Bezier.qBezier2(x1, y1, sx, sy, x2, y2, width_func, width_func_d);
      }
      else {
        if (a1 == STARTTYPE.THIN && a2 == ENDTYPE.OPEN) { // L2RD: fatten
          width_func = t => widfun(t, x1, y1, x2, y2, this.kMinWidthT) * this.kL2RDfatten;
          width_func_d = t => widfun_d(t, x1, y1, x2, y2, this.kMinWidthT) * this.kL2RDfatten;
        }
        else if (a1 == STARTTYPE.CONNECTING_V && a2 == ENDTYPE.OPEN) { //未使用。さんずい用 (実験)
          width_func = t => {return ((1-t)*0.628+Math.pow((1-t),30)*0.600+0.222)*this.kMinWidthT};
          width_func_d = t => {return (-0.628-30*Math.pow((1-t),29)*0.600)*this.kMinWidthT};
        }
        else if (a1 == STARTTYPE.THIN) {
          width_func = t => widfun(t, x1, y1, x2, y2, this.kMinWidthT);
          width_func_d = t => widfun_d(t, x1, y1, x2, y2, this.kMinWidthT);
        }
        else if (a2 == ENDTYPE.LEFT_SWEEP) {
          width_func = t => widfun(1 - t, x1, y1, x2, y2, this.kMinWidthT);
          width_func_d = t => -widfun_d(1 - t, x1, y1, x2, y2, this.kMinWidthT);
        }
        else {
          width_func = t => this.kMinWidthT;
          width_func_d = t => 0;
        }
        [bez1, bez2] = Bezier.qBezier(x1, y1, sx, sy, x2, y2, width_func, width_func_d);
      }
      if (a1 == 132 && x1 != sx) {
        let b1 = bezier_to_y(bez2[bez2.length - 1], y1);
        if (b1) { bez2[bez2.length - 1] = b1; }
        var temp = bez1[0].concat();//deep copy
        let b2 = bezier_to_y(temp.reverse(), y1);
        if (b2) { bez1[0] = b2.reverse(); }
      } else if (a1 == 22 && x1 != sx && y1 > y2) {
        let b1 = bezier_to_y(bez2[bez2.length - 1], y1);
        if (b1) { bez2[bez2.length - 1] = b1; }
        var temp = bez1[0].concat();//deep copy
        let b2 = bezier_to_y(temp.reverse(), y1 + 1);//??
        if (b2) { bez1[0] = b2.reverse(); }
      }
      var poly = Bezier.bez_to_poly(bez1);
      poly.concat(Bezier.bez_to_poly(bez2));
      cv.addPolygon(poly);

      const bez1e = bez1[bez1.length - 1][3];
      const bez1c2 = bez1[bez1.length - 1][2];

      const bez2s = bez2[0][0];
      const bez2c1 = bez2[0][1];
      const tan1 = [bez1e[0] - bez1c2[0], bez1e[1] - bez1c2[1]];
      const tan2 = [bez2s[0] - bez2c1[0], bez2s[1] - bez2c1[1]];
      return [tan1, tan2];
    }

    minchoDrawBezier(x1pre, y1pre, sx1, sy1, sx2, sy2, x2pre, y2pre, a1, a2, cv) {
      var delta;
      switch (a1) {
        case STARTTYPE.OPEN:
        case STARTTYPE.THIN:
          delta = -1 * this.kMinWidthY * 0.5;
          break;
        case STARTTYPE.UPPER_LEFT_CORNER:
          //case 32:
          delta = this.kMinWidthY;
          break;
        default:
          delta = 0;
          break;
      }
      let [x1, y1] = get_extended_dest(x1pre, y1pre, sx1, sy1, delta);

      switch (a2) {
        case ENDTYPE.STOP: // get shorten for tail's circle
          delta = -1 * this.kMinWidthT * 0.52;
          break;
        case ENDTYPE.TURN_LEFT:
          delta = -this.kMage;
          break;
        default:
          delta = 0;
          break;
      }
      let [x2, y2] = get_extended_dest(x2pre, y2pre, sx2, sy2, delta);

      var width_func;
      var width_func_d;
      let bez1, bez2;
      
      if (a1 == STARTTYPE.THIN && a2 == ENDTYPE.STOP) { //stop
              width_func = t => widfun_stop(t, x1, y1, x2, y2, this.kMinWidthT);
        width_func_d = t => widfun_stop_d(t, x1, y1, x2, y2, this.kMinWidthT);

        [bez1, bez2] = Bezier.cBezier(x1, y1, sx1, sy1, sx2, sy2, x2, y2, width_func, width_func_d);

        //width_func = t => widfun_fat(t, x1, y1, x2, y2, this.kMinWidthT);
        //width_func_d = t => widfun_fat_d(t, x1, y1, x2, y2, this.kMinWidthT);
        //[bez1, bez2] = Bezier.cBezier_slant(x1, y1, sx1, sy1, sx2, sy2, x2, y2, width_func, width_func_d);
      }
      else {
        if (a1 == STARTTYPE.THIN && a2 == ENDTYPE.OPEN) { // L2RD: fatten
          width_func = t => widfun(t, x1, y1, x2, y2, this.kMinWidthT) * this.kL2RDfatten;
          width_func_d = t => widfun_d(t, x1, y1, x2, y2, this.kMinWidthT) * this.kL2RDfatten;
        }
        else if (a1 == STARTTYPE.THIN) {
          width_func = t => widfun_fat(t, x1, y1, x2, y2, this.kMinWidthT);
          width_func_d = t => widfun_fat_d(t, x1, y1, x2, y2, this.kMinWidthT);
        }
        else if (a2 == ENDTYPE.LEFT_SWEEP) {
          width_func = t => widfun(1 - t, x1, y1, x2, y2, this.kMinWidthT);
          width_func_d = t => -widfun_d(1 - t, x1, y1, x2, y2, this.kMinWidthT);
        }
        else {
          width_func = t => this.kMinWidthT;
          width_func_d = t => 0;
        }
        [bez1, bez2] = Bezier.cBezier(x1, y1, sx1, sy1, sx2, sy2, x2, y2, width_func, width_func_d);
      }
      if (a1 == 132 && x1 != sx1) {
        let b1 = bezier_to_y(bez2[bez2.length - 1], y1);
        if (b1) { bez2[bez2.length - 1] = b1; }
        var temp = bez1[0];
        let b2 = bezier_to_y(temp.reverse(), y1);
        if (b2) { bez1[0] = b2.reverse(); }
      } else if (a1 == 22 && x1 > sx1) {
        let b1 = bezier_to_y(bez2[bez2.length - 1], y1);
        if (b1) { bez2[bez2.length - 1] = b1; }
        var temp = bez1[0];
        let b2 = bezier_to_y(temp.reverse(), y1 + 1);
        if (b2) { bez1[0] = b2.reverse(); }
      }
      var poly = Bezier.bez_to_poly(bez1);
      poly.concat(Bezier.bez_to_poly(bez2));
      cv.addPolygon(poly);

      const bez1e = bez1[bez1.length - 1][3];
      const bez1c2 = bez1[bez1.length - 1][2];

      const bez2s = bez2[0][0];
      const bez2c1 = bez2[0][1];
      const tan1 = [bez1e[0] - bez1c2[0], bez1e[1] - bez1c2[1]];
      const tan2 = [bez2s[0] - bez2c1[0], bez2s[1] - bez2c1[1]];
      return [tan1, tan2];
    }

    getStartOfVLine(x1, y1, x2, y2, a1, kMinWidthT, cv) {
      const dir = get_dir(x2 - x1, y2 - y1);
      var poly_start = new Polygon(2);
      if (dir.cos==0) {//vertical
        var left1, right1;
        switch (a1) {
          case 0:
            right1 = this.kMinWidthY / 2;
            left1 = -this.kMinWidthY / 2;
            break;
          case 12:
            right1 = this.kMinWidthY + kMinWidthT;
            left1 = this.kMinWidthY;
            break;
          case 32:
            right1 = this.kMinWidthY - 0.001;
            left1 = this.kMinWidthY - 0.001;
            break;
          case 1:
          case 6: //... no need
          case 22:
          default:
            right1 = 0;
            left1 = 0;
            break;
        }
        poly_start = this.getStartOfOffsetLine(x1, y1, dir, kMinWidthT, right1, left1);
        if (a1 == 22) { //box's right top corner
          cv.drawUpperRightCorner_straight_v(x1, y1, kMinWidthT, this.kMinWidthY, this.kWidth);
        }
        if (a1 == 0) { //beginning of the stroke
          cv.drawOpenBegin_straight(x1, y1, kMinWidthT, this.kMinWidthY, dir);
        }
      } else {
        const rad = Math.atan((y2 - y1) / (x2 - x1));
        const v = (x1 > x2) ? -1 : 1;
        if (a1 == 22) {
          if (dir.sin==0) {//error
            console.log("error: connecting_v at the end of the horizontal line");
            poly_start = this.getStartOfLine(x1, y1, dir, kMinWidthT);
          } else {
            poly_start.set(1, x1 + (kMinWidthT * v + 1) / Math.sin(rad), y1 + 1);//??
            poly_start.set(0, x1 - (kMinWidthT * v) / Math.sin(rad), y1);
          }
        } else if (a1 == 32) {
          if (dir.sin==0) {//error
            console.log("error: connecting_v at the end of the horizontal line");
            poly_start = this.getStartOfLine(x1, y1, dir, kMinWidthT);
          } else {
            poly_start.set(1, x1 + (kMinWidthT * v) / Math.sin(rad), y1);
            poly_start.set(0, x1 - (kMinWidthT * v) / Math.sin(rad), y1);
          }
        } else {
          var left1, right1;
          switch (a1) {
            case 0:
              right1 = this.kMinWidthY * 0.5;
              left1 = this.kMinWidthY * -0.5;
              break;
            case 12:
              right1 = this.kMinWidthY + kMinWidthT;
              left1 = this.kMinWidthY;
              break;
            case 1:
            case 6:
            default:
              right1 = 0;
              left1 = 0;
              break;
          }
          poly_start = this.getStartOfOffsetLine(x1, y1, dir, kMinWidthT, right1, left1);
        }
        if (a1 == 22) { //SHIKAKU MIGIUE UROKO NANAME DEMO MASSUGU MUKI
          cv.drawUpperRightCorner(x1, y1, kMinWidthT, this.kMinWidthY, this.kWidth);
        }
        if (a1 == 0) { //beginning of the stroke
          cv.drawOpenBegin_straight(x1, y1, kMinWidthT, this.kMinWidthY, dir);
        }
      }
      return poly_start;
    }

    getStartOfLine(x1, y1, dir, halfWidth) {
      //get polygon data for the start of line
      var poly = new Polygon(2);
      poly.set(1, x1 + dir.sin * halfWidth,
        y1 - dir.cos * halfWidth);
      poly.set(0, x1 - dir.sin * halfWidth,
        y1 + dir.sin * halfWidth);
      return poly;
    }

    getStartOfOffsetLine(x1, y1, dir, halfWidth, off_right1, off_left1) {
      //get polygon data for the start of line (with offset)
      var poly = new Polygon(2);
      poly.set(1, x1 + dir.sin * halfWidth - dir.cos * off_left1,
        y1 - dir.cos * halfWidth - dir.sin * off_left1);
      poly.set(0, x1 - dir.sin * halfWidth - dir.cos * off_right1,
        y1 + dir.cos * halfWidth - dir.sin * off_right1);
      return poly;
    }

    getEndOfLine(x1, y1, x2, y2, halfWidth) {
      //get polygon data for the end of line
      const dir = get_dir(x2 - x1, y2 - y1);
      var poly = new Polygon(2);
      poly.set(0, x2 + dir.sin * halfWidth,
        y2 - dir.cos * halfWidth);
      poly.set(1, x2 - dir.sin * halfWidth,
        y2 + dir.cos * halfWidth);
      return poly;
    }

    getEndOfOffsetLine(x1, y1, x2, y2, halfWidth, off_right2, off_left2) {
      //get polygon data for the end of line (with offset)
      const dir = get_dir(x2 - x1, y2 - y1);
      var poly = new Polygon(2);
      poly.set(0, x2 + dir.sin * halfWidth + off_left2 * dir.cos,
        y2 - dir.cos * halfWidth + off_left2 * dir.sin);
      poly.set(1, x2 - dir.sin * halfWidth + off_right2 * dir.cos,
        y2 + dir.cos * halfWidth + off_right2 * dir.sin);
      return poly;
    }

    //functions for adjustment

    adjustTateParam(stroke, others) { // strokes
      //(STROKETYPE.STRAIGHT || STROKETYPE.BENDING || STROKETYPE.VCURVE)
      if (stroke[3] != stroke[5]) return 0;
      var res = 0;
      for (let other of others) {
        if ((other[0] == 1 || other[0] == 3 || other[0] == 7) && other[3] == other[5] &&
          !(stroke[4] + 1 > other[6] || stroke[6] - 1 < other[4]) &&
          Math.abs(stroke[3] - other[3]) < this.kMinWidthT * this.kAdjustTateStep) {
          res += (this.kAdjustTateStep - Math.floor(Math.abs(stroke[3] - other[3]) / this.kMinWidthT));
        }
      }
      res = Math.min(res, this.kAdjustTateStep);
      return res;//a2 += res * 1000
    }

    adjustUrokoParam(stroke, others) { // strokes
      //STROKETYPE.STRAIGHT && ENDTYPE.OPEN
      for (var k = 0; k < this.kAdjustUrokoLengthStep; k++) {
        var tx, ty, tlen;
        if (stroke[4] == stroke[6]) { // YOKO
          tx = stroke[5] - this.kAdjustUrokoLine[k];
          ty = stroke[6] - 0.5;
          tlen = stroke[5] - stroke[3];
        } else {
          var rad = Math.atan((stroke[6] - stroke[4]) / (stroke[5] - stroke[3]));
          tx = stroke[5] - this.kAdjustUrokoLine[k] * Math.cos(rad) - 0.5 * Math.sin(rad);
          ty = stroke[6] - this.kAdjustUrokoLine[k] * Math.sin(rad) - 0.5 * Math.cos(rad);
          tlen = Math.sqrt((stroke[6] - stroke[4]) * (stroke[6] - stroke[4]) +
            (stroke[5] - stroke[3]) * (stroke[5] - stroke[3]));
        }
        if (tlen < this.kAdjustUrokoLength[k] ||
          isCrossWithOthers(others, -1, tx, ty, stroke[5], stroke[6])
        ) {
          return (this.kAdjustUrokoLengthStep - k);
        }
      }
      return 0;//a3 += res * 100;
    }

    adjustUroko2Param(stroke, others) { // strokes
      //STROKETYPE.STRAIGHT && ENDTYPE.OPEN && y1==y2
      var pressure = 0;
      for (let other of others) {
        if (
          (other[0] == 1 && other[4] == other[6] &&
            !(stroke[3] + 1 > other[5] || stroke[5] - 1 < other[3]) &&
            Math.abs(stroke[4] - other[4]) < this.kAdjustUroko2Length) ||
          (other[0] == 3 && other[6] == other[8] &&
            !(stroke[3] + 1 > other[7] || stroke[5] - 1 < other[5]) &&
            Math.abs(stroke[4] - other[6]) < this.kAdjustUroko2Length)
        ) {
          pressure += Math.pow(this.kAdjustUroko2Length - Math.abs(stroke[4] - other[6]), 1.1);
        }
      }
      var result = Math.min(Math.floor(pressure / this.kAdjustUroko2Length), this.kAdjustUroko2Step);
      return result;//a3 += res * 100;
    }

    adjustHaneParam(epx, epy, others) { // adjust "Hane" (short line turning to the left)
      //endPointX, endPointY
      //(STROKETYPE.STRAIGHT || STROKETYPE.CURVE || STROKETYPE.BEZIER) && ENDTYPE.TURN_LEFT
      var res = 0;
      var nearest = Infinity; // the nearest point to the short line
      if (epx + 18 < 100) {
        nearest = epx + 18;
      }
      for (let other of others) {
        if (other[0] == STROKETYPE.STRAIGHT && other[3] == other[5] && other[3] < epx && other[4] <= epy && other[6] >= epy) {
          if (epx - other[3] < 100) {
            nearest = Math.min(nearest, epx - other[3]);
          }
        }
      }
      if (nearest != Infinity) {
        res = 7 - Math.floor(nearest / 15);
      }
      return res;//a3 += res * 100;
    }

    adjustMageParam(stroke, others) {
      //STROKETYPE.BENDING
      //applied only if y2=y3
      if (stroke[6] != stroke[8]) return 0;
      var res = 0;
      for (let other of others) {
        if (
          (other[0] == 1 && other[4] == other[6] &&
            !(stroke[5] + 1 > other[5] || stroke[7] - 1 < other[3]) &&
            Math.abs(stroke[6] - other[4]) < this.kMinWidthT * this.kAdjustMageStep) ||
          (other[0] == 3 && other[6] == other[8] &&
            !(stroke[5] + 1 > other[7] || stroke[7] - 1 < other[5]) &&
            Math.abs(stroke[6] - other[6]) < this.kMinWidthT * this.kAdjustMageStep)
        ) {
          res += this.kAdjustMageStep - Math.floor(Math.abs(stroke[6] - other[6]) / this.kMinWidthT);
        }
      }
      res = Math.min(res, this.kAdjustMageStep);
      return res;//a3 += res * 1000;
    }

    adjustKirikuchiParam(stroke, others) { // connecting to other strokes.
      //STROKETYPE.CURVE, STARTTYPE.CONNECTING_V
      if (stroke[3] > stroke[5] &&
        stroke[4] < stroke[6]) {
        for (let other of others) {
          if (other[0] == 1 &&
            other[3] < stroke[3] && other[5] > stroke[3] &&
            other[4] == stroke[4] && other[4] == other[6]) {
            return true;
          }
        }
      }
      return false;
      //if (res) a2 += 100;
    }
    
    adjustKakatoParam(stroke, others) {
      //if (STROKETYPE.STRAIGHT && (LOWER_LEFT_CORNER || LOWER_RIGHT_CORNER))
      for (var k = 0; k < this.kAdjustKakatoStep; k++) {
        if (isCrossBoxWithOthers(others, -1,
          stroke[5] - this.kAdjustKakatoRangeX / 2,
          stroke[6] + this.kAdjustKakatoRangeY[k],
          stroke[5] + this.kAdjustKakatoRangeX / 2,
          stroke[6] + this.kAdjustKakatoRangeY[k + 1])
          | stroke[6] + this.kAdjustKakatoRangeY[k + 1] > 200 // adjust for baseline
          | stroke[6] - stroke[4] < this.kAdjustKakatoRangeY[k + 1] // for thin box
        ) {
          return (3 - k);
        }
      }
      return 0;
      //a3 += res * 100;
    }
  }

  const FONTTYPE = {
    MINCHO: 0,
    GOTHIC: 1,
  };

  class Kage {
    constructor(type, size){
      this.kBuhin = new Buhin();
      this.setFont(type,size);
      this.kRate = 100;
    }
    setFont(type, size){
      switch(type){
       case FONTTYPE.GOTHIC:{
          this.kFont = new Gothic(size);
          break;
        }
        case FONTTYPE.MINCHO:{
          this.kFont = new Mincho(size);
          break;
        }
        default:{
          this.kFont = new Mincho(size);
          break;
        }
      }
    }
    makeGlyph(polygons, buhin) { // The word "buhin" means "component".  This method converts buhin (KAGE data format) to polygons (path data).  The variable buhin may represent a component of kanji or a kanji itself.
      var glyphData = this.kBuhin.search(buhin);
      this.makeGlyph2(polygons, glyphData);
    }
    makeGlyph2(polygons, data) {
        var kageStrokes = this.getStrokes(data);
        polygons.concat(this.kFont.getPolygons(kageStrokes));
    }
    makeGlyph3(data) { // void
      var kageStrokes = this.getStrokes(data);
      return this.kFont.getPolygons(kageStrokes);
    }
    getStrokes(glyphData) { // strokes array
      var strokes = new Array();
      var textData = glyphData.split("$");
      for (var i = 0; i < textData.length; i++) {
        var columns = textData[i].split(":");
        if (Math.floor(columns[0]) != STROKETYPE.REFERENCE) {
          strokes.push([
            Math.floor(columns[0]),
            Math.floor(columns[1]),
            Math.floor(columns[2]),
            Math.floor(columns[3]),
            Math.floor(columns[4]),
            Math.floor(columns[5]),
            Math.floor(columns[6]),
            Math.floor(columns[7]),
            Math.floor(columns[8]),
            Math.floor(columns[9]),
            Math.floor(columns[10])
          ]);

        } else {
          var buhin = this.kBuhin.search(columns[7]);
          if (buhin != "") {
            strokes = strokes.concat(this.getStrokesOfBuhin(buhin,
              Math.floor(columns[3]),
              Math.floor(columns[4]),
              Math.floor(columns[5]),
              Math.floor(columns[6]),
              Math.floor(columns[1]),
              Math.floor(columns[2]),
              Math.floor(columns[9]),
              Math.floor(columns[10]))
            );
          }
        }
      }
      return strokes;
    }

    getStrokesOfBuhin(buhin, x1, y1, x2, y2, sx, sy, sx2, sy2) {
      var temp = this.getStrokes(buhin);
      var result = new Array();
      var box = getBoundingBox(temp);
      if (sx != 0 || sy != 0) {
        if (sx > 100) {
          sx -= 200;
        } else {
          sx2 = 0;
          sy2 = 0;
        }
      }
      for (var i = 0; i < temp.length; i++) {
        if (sx != 0 || sy != 0) {
          temp[i][3] = stretch(sx, sx2, temp[i][3], box.minX, box.maxX);
          temp[i][4] = stretch(sy, sy2, temp[i][4], box.minY, box.maxY);
          temp[i][5] = stretch(sx, sx2, temp[i][5], box.minX, box.maxX);
          temp[i][6] = stretch(sy, sy2, temp[i][6], box.minY, box.maxY);
          if (temp[i][0] != STROKETYPE.REFERENCE) {
            temp[i][7] = stretch(sx, sx2, temp[i][7], box.minX, box.maxX);
            temp[i][8] = stretch(sy, sy2, temp[i][8], box.minY, box.maxY);
            temp[i][9] = stretch(sx, sx2, temp[i][9], box.minX, box.maxX);
            temp[i][10] = stretch(sy, sy2, temp[i][10], box.minY, box.maxY);
          }
        }
        result.push([temp[i][0],
        temp[i][1],
        temp[i][2],
        x1 + temp[i][3] * (x2 - x1) / 200,
        y1 + temp[i][4] * (y2 - y1) / 200,
        x1 + temp[i][5] * (x2 - x1) / 200,
        y1 + temp[i][6] * (y2 - y1) / 200,
        x1 + temp[i][7] * (x2 - x1) / 200,
        y1 + temp[i][8] * (y2 - y1) / 200,
        x1 + temp[i][9] * (x2 - x1) / 200,
        y1 + temp[i][10] * (y2 - y1) / 200]);
        
      }
      return result;
    }
   
  }

  exports.FONTTYPE = FONTTYPE;
  exports.Kage = Kage;
  exports.Polygons = Polygons;

  Object.defineProperty(exports, '__esModule', { value: true });

  return exports;

})({});
