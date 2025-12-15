import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  base: './',
  build: {
    sourcemap: true,
  },
  worker: {
    format: 'es',
  },

  // Vite serves .gz files with Content-Encoding: gzip
  // which the browser auto-decompresses, making JBrowse not work.
  // The only work-around is to serve it separately without vite
  // in development mode.
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
