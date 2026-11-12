"use strict";

const svg = document.querySelector("#chart svg");
const geometry = JSON.parse(document.querySelector("#sector-geometry").textContent);
const selection = document.querySelector("#selection");
const clear = document.querySelector("#clear");
const scaleLabel = document.querySelector("#scale");
const view = svg.viewBox.baseVal;
const initial = { x: view.x, y: view.y, width: view.width, height: view.height };
let box = { ...initial };
let drag = null;
let pendingFrame = null;
let bounds = null;
let selected = null;

// Reuse one overlay. Point selection does not rebuild the chart or its ribbons.
const highlight = document.createElementNS("http://www.w3.org/2000/svg", "path");
highlight.id = "sector-highlight";
highlight.setAttribute("aria-hidden", "true");
svg.appendChild(highlight);

function schedulePaint() {
  if (pendingFrame !== null) return;
  pendingFrame = requestAnimationFrame(() => {
    pendingFrame = null;
    bounds = null;
    svg.setAttribute("viewBox", [box.x, box.y, box.width, box.height].join(" "));
    scaleLabel.textContent = Math.round(initial.width / box.width * 100) + "%";
  });
}

function viewport() {
  if (!bounds) bounds = svg.getBoundingClientRect();
  const scale = Math.min(bounds.width / box.width, bounds.height / box.height);
  if (!Number.isFinite(scale) || scale <= 0) return null;
  return {
    scale,
    left: bounds.left + (bounds.width - box.width * scale) / 2,
    top: bounds.top + (bounds.height - box.height * scale) / 2
  };
}

function pointerPoint(event) {
  const v = viewport();
  return v ? {
    x: box.x + (event.clientX - v.left) / v.scale,
    y: box.y + (event.clientY - v.top) / v.scale
  } : null;
}

function zoom(factor, anchor) {
  const ratio = Math.max(0.02, Math.min(2, box.width * factor / initial.width));
  const nextWidth = initial.width * ratio;
  const nextFactor = nextWidth / box.width;
  const point = anchor || { x: box.x + box.width / 2, y: box.y + box.height / 2 };
  box = {
    x: point.x + (box.x - point.x) * nextFactor,
    y: point.y + (box.y - point.y) * nextFactor,
    width: nextWidth,
    height: box.height * nextFactor
  };
  schedulePaint();
}

document.querySelector("#plus").addEventListener("click", () => zoom(0.8));
document.querySelector("#minus").addEventListener("click", () => zoom(1.25));
document.querySelector("#fit").addEventListener("click", () => {
  box = { ...initial };
  schedulePaint();
});
svg.addEventListener("wheel", event => {
  if (!event.deltaY) return;
  event.preventDefault();
  zoom(event.deltaY < 0 ? 0.9 : 1.1, pointerPoint(event));
}, { passive: false });

svg.addEventListener("pointerdown", event => {
  if (event.button !== 0 || drag) return;
  const v = viewport();
  if (!v) return;
  drag = { id: event.pointerId, x: event.clientX, y: event.clientY };
  svg.setPointerCapture(event.pointerId);
});
svg.addEventListener("pointermove", event => {
  if (!drag || drag.id !== event.pointerId) return;
  const v = viewport();
  if (!v) return;
  box.x -= (event.clientX - drag.x) / v.scale;
  box.y -= (event.clientY - drag.y) / v.scale;
  drag.x = event.clientX;
  drag.y = event.clientY;
  schedulePaint();
});
for (const type of ["pointerup", "pointercancel", "lostpointercapture"]) {
  svg.addEventListener(type, event => {
    if (drag && drag.id === event.pointerId) drag = null;
  });
}
function invalidateBounds() { bounds = null; }
window.addEventListener("resize", invalidateBounds);
window.addEventListener("scroll", invalidateBounds, true);

function clearSelection() {
  if (selected) selected.setAttribute("aria-pressed", "false");
  selected = null;
  highlight.removeAttribute("d");
  selection.textContent = "";
  clear.hidden = true;
}

document.querySelector("aside").addEventListener("click", event => {
  const entry = event.target.closest(".sector-entry");
  if (!entry) return;
  const alreadySelected = entry === selected;
  clearSelection();
  if (alreadySelected) return;
  const sector = geometry[Number(entry.dataset.sector) - 1];
  if (!sector) return;
  selected = entry;
  selected.setAttribute("aria-pressed", "true");
  highlight.setAttribute("d", sector.path);
  selection.textContent = entry.dataset.name + " (" + entry.dataset.rank + ")";
  clear.hidden = false;
  const [x, y] = sector.center;
  // Reveal a selected sector if the current zoom has panned away from it.
  if (x < box.x || x > box.x + box.width || y < box.y || y > box.y + box.height) {
    box.x = x - box.width / 2;
    box.y = y - box.height / 2;
    schedulePaint();
  }
});
clear.addEventListener("click", clearSelection);
document.addEventListener("keydown", event => {
  if (event.key === "Escape") clearSelection();
});

const search = document.querySelector("#search");
const rows = Array.from(document.querySelectorAll(".sector-row"), row => ({
  element: row, text: row.textContent.toLowerCase()
}));
let searchTimer;
search.addEventListener("input", () => {
  clearTimeout(searchTimer);
  searchTimer = setTimeout(() => {
    const query = search.value.trim().toLowerCase();
    for (const row of rows) row.element.hidden = !row.text.includes(query);
    if (query) document.querySelectorAll("aside details").forEach(details => { details.open = true; });
  }, 100);
});
schedulePaint();
