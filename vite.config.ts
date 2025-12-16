import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  base: './', // ./ means assets/data are relative to index.html
  worker: {
    format: 'es',
  },

  // keep index.html in ./src/ (Vite sucks)
  root: 'src',
  publicDir: '../public',
  build: {
    outDir: '../dist',
    emptyOutDir: true,
    sourcemap: true,
  },

  // Vite serves .gz files with "Content-Encoding: gzip" which
  // the browser auto-decompresses, making JBrowse not work!!!
  // The only work-around is to serve it separately without vite
  // in dev mode (static-build is unaffected)
  server: {
    proxy: {
      '/data': {
        target: 'http://localhost:5174',
        changeOrigin: true,
        rewrite: p => p.replace(/^\/data/, ''),
      },
    },
  },
});
