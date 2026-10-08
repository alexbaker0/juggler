// Scene-file parser and camera setup for the Juggler ray tracer.
// Port of load_scene() from juggler_orig.c (the PC port of Eric Graham's
// RT1.C/RT2.C); the ray tracing itself is in rt.wasm (src/rt.c).
// Used by worker.js and by Node tests.

// Parse a robot.dat-style scene. Throws Error with a readable message on bad input.
function parseScene(text) {
  let p = 0;
  const fail = what => { throw new Error(`${what} near "${text.slice(p, p + 20).replace(/\s+/g, " ")}"`); };
  const skipws = () => { while (p < text.length && " \t\n\r".includes(text[p])) p++; };
  const peek = () => { skipws(); return p < text.length ? text[p] : ""; };
  const expect = c => { if (peek() !== c) fail(`expected '${c}'`); p++; };
  // A number, skipping whitespace and the , and : field separators (like p_num).
  const num = () => {
    skipws();
    while (p < text.length && (text[p] === "," || text[p] === ":")) { p++; skipws(); }
    const m = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?/.exec(text.slice(p));
    if (!m) fail("expected a number");
    p += m[0].length;
    return parseFloat(m[0]);
  };
  const isNumStart = c => /[0-9+\-.]/.test(c);
  const point = () => { expect("("); const d = [num(), num(), num()]; expect(")"); return d; };
  const color = () => { expect("<"); const d = [num(), num(), num()]; expect(">"); return d; };

  const eye = point();
  expect("["); const alt = num(), az = num(); expect("]");
  const lens = num();

  // Objects: <r,g,b> type (x,y,z):rad [n (x,y,z):rad]... ;  -- each "n (..):r"
  // sweeps n spheres from the previous point/radius (a limb). ';;' ends the list.
  const spheres = [];
  while (peek() === "<") {
    const col = color(), type = Math.trunc(num());
    let cur = point(); expect(":"); let curR = num();
    spheres.push({ pos: cur, rad: curR, col, type });
    while (isNumStart(peek())) {
      const count = Math.trunc(num());
      const nxt = point(); expect(":"); const nxtR = num();
      for (let i = 1; i <= count; i++) {
        const t = i / count;
        spheres.push({ pos: [0, 1, 2].map(k => cur[k] + (nxt[k] - cur[k]) * t), rad: curR + (nxtR - curR) * t, col, type });
      }
      cur = nxt; curR = nxtR;
    }
    expect(";");
  }
  while (peek() === ";") p++;

  const lamps = [];
  const nl = Math.trunc(num());
  for (let i = 0; i < nl; i++) {
    const pos = point(); expect(":"); const rad = num();
    lamps.push({ pos, rad, col: color() });
  }
  // Five world colours: tile0, tile1, ambient illum, sky zenith, sky horizon.
  const [tile0, tile1, illum, skyZen, skyHor] = [color(), color(), color(), color(), color()];
  if (spheres.length > 2048) throw new Error(`too many spheres (${spheres.length}, max 2048)`);
  return { eye, alt, az, lens, spheres, lamps, tile0, tile1, illum, skyZen, skyHor };
}

// Load a parsed scene into a rt.wasm instance's exports, for an nx x ny image.
// altDeg/azDeg override the scene's view angles (for the camera controls).
function setupScene(ex, s, nx, ny, altDeg = s.alt, azDeg = s.az) {
  ex.clear_scene();
  for (const o of s.spheres) ex.add_sphere(...o.pos, o.rad, ...o.col, o.type);
  for (const l of s.lamps) ex.add_lamp(...l.pos, l.rad, ...l.col);
  ex.set_world(...s.tile0, ...s.tile1, ...s.illum, ...s.skyZen, ...s.skyHor);
  // Camera exactly as RT2.C: degrees -> radians with the same pi constant.
  const alt = altDeg * (3.14159265358979 / 180.0), az = azDeg * (3.14159265358979 / 180.0);
  ex.set_camera(...s.eye,
    Math.cos(az) * Math.cos(alt), Math.sin(az) * Math.cos(alt), Math.sin(alt),
    Math.sin(az), -Math.cos(az), 0.0,
    -Math.cos(az) * Math.sin(alt), -Math.sin(az) * Math.sin(alt), Math.cos(alt),
    0.028 * s.lens, nx, ny);
}

if (typeof module !== "undefined") module.exports = { parseScene, setupScene };
