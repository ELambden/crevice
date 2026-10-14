/* CREVICE interactive scene loader for 3Dmol.js (MIT licence, part of CREVICE).
 * Builds the layers of a crevice.interactive_view scene into a 3Dmol viewer
 * and, when a container is given, one checkbox per layer. No text is drawn
 * into the 3D canvas. */
(function (root) {
  "use strict";

  function bytesFromBase64(text) {
    var raw = atob(text), out = new Uint8Array(raw.length);
    for (var i = 0; i < raw.length; i++) out[i] = raw.charCodeAt(i);
    return out;
  }

  function decodeArray(spec) {
    var bytes = bytesFromBase64(spec.data);
    var done = function (buffer) {
      return spec.dtype === "uint32" ? new Uint32Array(buffer) : new Float32Array(buffer);
    };
    if (spec.encoding === "gzip+base64") {
      var stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
      return new Response(stream).arrayBuffer().then(done);
    }
    return Promise.resolve(done(bytes.buffer));
  }

  function points(array) {
    var out = [];
    for (var i = 0; i + 2 < array.length; i += 3) out.push({x: array[i], y: array[i + 1], z: array[i + 2]});
    return out;
  }

  function meshSpec(layer) {
    return Promise.all([decodeArray(layer.vertices), decodeArray(layer.faces), decodeArray(layer.normals)])
      .then(function (arrays) {
        return {vertexArr: points(arrays[0]), faceArr: Array.from(arrays[1]), normalArr: points(arrays[2]),
                color: layer.color, opacity: layer.opacity};
      });
  }

  function build(viewer, controls, scene) {
    var model = null, state = {}, shapes = {}, specs = {};
    var layers = scene.layers;
    if (scene.structure) model = viewer.addModel(scene.structure.data, scene.structure.format);

    function restyle() {
      if (!model) return;
      viewer.setStyle({model: model}, {});
      layers.forEach(function (layer) {
        if (!state[layer.name]) return;
        if (layer.kind === "cartoon") {
          viewer.addStyle({model: model}, {cartoon: {color: layer.color, opacity: layer.opacity}});
        } else if (layer.kind === "sticks") {
          var map = {C: layer.color, N: "#3050f8", O: "#ff0d0d", S: "#e6c200", P: "#ff8000", H: "#ffffff"};
          layer.residues.forEach(function (group) {
            var sel = {model: model, resi: group.resi};
            if (group.chain) sel.chain = group.chain;
            viewer.addStyle(sel, {stick: {radius: layer.radius, colorscheme: {prop: "elem", map: map}}});
          });
        }
      });
    }

    function drawShape(layer) {
      if (layer.kind === "mesh") {
        shapes[layer.name] = [viewer.addCustom(specs[layer.name])];
      } else if (layer.kind === "line") {
        shapes[layer.name] = [];
        layer.paths.forEach(function (path) {
          for (var i = 0; i + 1 < path.length; i++) {
            shapes[layer.name].push(viewer.addCylinder({
              start: {x: path[i][0], y: path[i][1], z: path[i][2]},
              end: {x: path[i + 1][0], y: path[i + 1][1], z: path[i + 1][2]},
              radius: layer.radius, color: layer.color, fromCap: 1, toCap: 1}));
          }
        });
      } else if (layer.kind === "spheres") {
        shapes[layer.name] = layer.points.map(function (p) {
          return viewer.addSphere({center: {x: p[0], y: p[1], z: p[2]}, radius: layer.radius,
                                   color: layer.color, opacity: layer.opacity});
        });
      }
    }

    function setVisible(name, visible) {
      var layer = layers.find(function (l) { return l.name === name; });
      if (!layer) return;
      state[name] = !!visible;
      if (layer.kind === "cartoon" || layer.kind === "sticks") {
        restyle();
      } else {
        (shapes[name] || []).forEach(function (s) { viewer.removeShape(s); });
        shapes[name] = [];
        if (visible) drawShape(layer);
      }
      viewer.render();
    }

    var pending = layers.filter(function (l) { return l.kind === "mesh"; }).map(function (layer) {
      return meshSpec(layer).then(function (spec) { specs[layer.name] = spec; });
    });

    return Promise.all(pending).then(function () {
      layers.forEach(function (layer) {
        state[layer.name] = !!layer.visible;
        if (layer.visible && layer.kind !== "cartoon" && layer.kind !== "sticks") drawShape(layer);
      });
      restyle();
      viewer.zoomTo();
      if (scene.camera && scene.camera.quaternion) {
        var v = viewer.getView();
        var c = scene.camera.center;
        v[0] = -c[0]; v[1] = -c[1]; v[2] = -c[2];
        v[4] = scene.camera.quaternion[0]; v[5] = scene.camera.quaternion[1];
        v[6] = scene.camera.quaternion[2]; v[7] = scene.camera.quaternion[3];
        viewer.setView(v);
        if (scene.camera.zoom) viewer.zoom(scene.camera.zoom);
      }
      var home = viewer.getView();
      viewer.render();
      if (controls) buildControls(controls, layers, state, setVisible, viewer, home);
      return {viewer: viewer, setVisible: setVisible, state: state};
    });
  }

  function buildControls(controls, layers, state, setVisible, viewer, home) {
    controls.innerHTML = "";
    layers.forEach(function (layer) {
      var row = document.createElement("label");
      row.className = "crevice-layer";
      var box = document.createElement("input");
      box.type = "checkbox";
      box.checked = state[layer.name];
      box.dataset.layer = layer.name;
      box.addEventListener("change", function () { setVisible(layer.name, box.checked); });
      var swatch = document.createElement("span");
      swatch.className = "crevice-swatch";
      swatch.style.background = layer.color;
      var text = document.createElement("span");
      text.textContent = layer.label;
      if (layer.description) row.title = layer.description;
      row.appendChild(box); row.appendChild(swatch); row.appendChild(text);
      controls.appendChild(row);
    });
    var buttons = document.createElement("div");
    buttons.className = "crevice-buttons";
    var reset = document.createElement("button");
    reset.type = "button"; reset.textContent = "Reset view";
    reset.addEventListener("click", function () { viewer.setView(home); viewer.render(); });
    var png = document.createElement("button");
    png.type = "button"; png.textContent = "Save PNG";
    png.addEventListener("click", function () {
      var a = document.createElement("a");
      a.href = viewer.pngURI(); a.download = "crevice_view.png"; a.click();
    });
    buttons.appendChild(reset); buttons.appendChild(png);
    controls.appendChild(buttons);
  }

  root.CreviceView = {build: build, decodeArray: decodeArray};
})(typeof window !== "undefined" ? window : this);
