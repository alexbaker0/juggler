# The Juggler

Eric Graham's 1986 Amiga ray tracer, rendering his original Juggler scene live in the browser.

**Try it:** https://alexbaker0.github.io/juggler/

- `src/rt.c`: Graham's RT1.C engine (Copyright 1987 Eric Graham; "Permission is granted to copy and modify this file, provided that this notice is retained"). This is the PC port with ANSI prototypes, a fix to `reflect()` (the scanned original divided by the normal), and a mirror-recursion guard, plus a small row-rendering API. Compiled to `rt.wasm`.
- `scene.js`: the `robot.dat` parser and RT2.C camera/lighting setup, ported to JS.
- `worker.js`: rows are split across Web Workers, one per CPU core, and drawn as they finish.
- At 1024×768 with 3×3 anti-aliasing the output is pixel-identical to the native Windows build.

Build: `clang --target=wasm32 -O3 -mbulk-memory -nostdlib -Wl,--no-entry -o rt.wasm src/rt.c`
