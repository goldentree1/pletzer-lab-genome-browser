import Plugin from '@jbrowse/core/Plugin';
import PluginManager from '@jbrowse/core/PluginManager';

export default class PluginTest extends Plugin {
  name = 'PluginTest';

  install(pluginManager: PluginManager) {
    // sanity check: proves plugin is loaded
    console.log('[PluginTest] installed');

    // example: remove track menu for GFF3 FeatureTracks
    pluginManager.addToExtensionPoint('TrackMenuItems', (items, ctx) => {
      const model = ctx?.model;

      if (
        /** @ts-expect-error fuck typescript */
        model?.configuration?.trackType === 'FeatureTrack' &&
        /** @ts-expect-error fuck typescript #2 */
        model?.adapterConfig?.type === 'Gff3TabixAdapter'
      ) {
        return []; // nukes the ⋮ menu
      }

      pluginManager.addToExtensionPoint('FeatureMenuItems', (items, ctx) => {
        /** @ts-expect-error fuck typescript */
        if (ctx?.track?.adapterConfig?.type === 'Gff3TabixAdapter') {
          return [];
        }
        return items;
      });

      return items;
    });
  }
}
