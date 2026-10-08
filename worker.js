// One render worker. Gets the compiled rt.wasm module, the scene text, image
// size, anti-aliasing and camera angles, plus which rows to render; posts each
// finished RGBA row back as it goes.
importScripts("scene.js");

onmessage = async ({ data }) => {
  const { module, sceneText, nx, ny, aa, alt, az, eye, rows } = data;
  const ex = (await WebAssembly.instantiate(module, {})).exports;
  const s = parseScene(sceneText);
  if (eye) s.eye = eye;
  setupScene(ex, s, nx, ny, alt, az);
  const row = new Uint8Array(ex.memory.buffer, ex.row_buffer(), nx * 4);
  for (const j of rows) {
    ex.render_row(j, aa);
    const copy = row.slice();
    postMessage({ j, pixels: copy }, [copy.buffer]);
  }
  postMessage({ done: true });
};
