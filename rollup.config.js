export default
  [
    {
      input: 'index.js',
      output: {
        file: 'dist/kage_bundle.js',
        format: 'iife',
        name: 'kage_export'
      }
    },
    {
      input: 'index.js',
      output: {
        file: 'dist/kage_bundle_cjs.js',
        format: 'cjs'
      }
    }
  ];