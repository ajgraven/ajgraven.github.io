// =============================================================================
// sw.js -- Service worker for offline / PWA support (P3.3).
//
// On install: pre-cache the static asset list from asset-manifest.js.
// On fetch:
//   * same-origin → cache-first (offline-safe; latest version after a load).
//   * CDN (math.js, KaTeX) → network-first with cache fallback.
// On activate: drop any stale caches whose name doesn't match CACHE_NAME.
//
// Cache versioning: every release that ships changes to app/ should bump
// CACHE_VERSION below. Old caches are dropped automatically on activate,
// forcing a re-fetch on next load.
// =============================================================================

'use strict';

// Import the shared asset manifest. importScripts is the only way to load
// helpers from a classic service-worker scope.
importScripts('./asset-manifest.js');

// Single canonical version lives in asset-manifest.js (imported above), so a
// release bumps it in exactly one place. Because that file is importScripts'd
// here, bumping the constant is part of this worker's byte-update comparison,
// which reliably triggers a service-worker update → old caches dropped on
// activate (see :48-53).
const CACHE_VERSION = self.QD_ASSET_MANIFEST.CACHE_VERSION;
const CACHE_NAME    = 'qd-solver-' + CACHE_VERSION;

const { ALL_ASSETS, CDN_ASSETS } = self.QD_ASSET_MANIFEST;

// ---------------------------------------------------------------------------
// Install: pre-cache the static asset list. CDN assets are NOT pre-cached
// (the page may load offline and we don't want install to fail) — they
// land in the cache on first online visit via the fetch handler.
// ---------------------------------------------------------------------------
self.addEventListener('install', (event) => {
  event.waitUntil((async () => {
    const cache = await caches.open(CACHE_NAME);
    // Use individual addAll so one missing file doesn't fail the whole install.
    await Promise.all(ALL_ASSETS.map(async (url) => {
      try { await cache.add(url); }
      catch (e) { console.warn('[sw] skip pre-cache for ' + url + ': ' + (e && e.message || e)); }
    }));
    // Take control of any open clients on next reload.
    self.skipWaiting();
  })());
});

// ---------------------------------------------------------------------------
// Activate: drop any older caches.
// ---------------------------------------------------------------------------
self.addEventListener('activate', (event) => {
  event.waitUntil((async () => {
    const names = await caches.keys();
    await Promise.all(names.map((name) => {
      if (name.startsWith('qd-solver-') && name !== CACHE_NAME) {
        return caches.delete(name);
      }
    }));
    await self.clients.claim();
  })());
});

// ---------------------------------------------------------------------------
// Fetch: cache-first for same-origin; network-first for CDN.
// ---------------------------------------------------------------------------
self.addEventListener('fetch', (event) => {
  const req = event.request;
  // Only handle GET; PUTs / POSTs / etc. fall through to default behavior.
  if (req.method !== 'GET') return;

  const url = new URL(req.url);
  const isSameOrigin = (url.origin === self.location.origin);
  const isCDN = CDN_ASSETS.some((a) => req.url === a);

  if (isSameOrigin) {
    event.respondWith(cacheFirst(req));
  } else if (isCDN) {
    event.respondWith(networkFirst(req));
  }
  // Other cross-origin requests: untouched.
});

async function cacheFirst(req) {
  const cache = await caches.open(CACHE_NAME);
  // ignoreSearch: the page loads scripts as `foo.js?v=<CACHE_VERSION>` but the
  // install step pre-caches the bare paths. Matching with ignoreSearch lets the
  // precache actually HIT (otherwise every request was a cache miss → needless
  // re-fetch). This is safe because freshness across deploys is enforced by the
  // versioned CACHE_NAME (a new release → new cache → activate drops the old),
  // NOT by the `?v=` query — within one release every `?v=` value is identical.
  const cached = await cache.match(req, { ignoreSearch: true });
  if (cached) return cached;
  try {
    const fresh = await fetch(req);
    // Only cache successful responses; opaque (CORS-blocked) responses are skipped.
    if (fresh && fresh.ok) cache.put(req, fresh.clone());
    return fresh;
  } catch (e) {
    // Offline + not in cache → return whatever we have (or 504).
    return cached || new Response('Offline and not cached', {
      status: 504, statusText: 'Gateway Timeout',
    });
  }
}

async function networkFirst(req) {
  const cache = await caches.open(CACHE_NAME);
  try {
    const fresh = await fetch(req);
    if (fresh && (fresh.ok || fresh.type === 'opaque')) {
      // Cache CDN responses (including opaque) so offline reloads work.
      cache.put(req, fresh.clone());
    }
    return fresh;
  } catch (e) {
    const cached = await cache.match(req);
    if (cached) return cached;
    return new Response('Offline and not cached: ' + req.url, {
      status: 504, statusText: 'Gateway Timeout',
    });
  }
}
