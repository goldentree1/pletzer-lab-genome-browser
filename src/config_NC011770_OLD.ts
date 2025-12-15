const conf = {
  assemblies: [
    {
      name: 'assembly1',
      sequence: {
        type: 'ReferenceSequenceTrack',
        trackId: 'FastaReferenceSequenceTrack',
        adapter: {
          type: 'BgzipFastaAdapter',
          fastaLocation: {
            uri: '/data/NC011770/LESB58_ASM2664v1.fasta.gz',
            locationType: 'UriLocation',
          },
          faiLocation: {
            uri: '/data/NC011770/LESB58_ASM2664v1.fasta.gz.fai',
            locationType: 'UriLocation',
          },
          gziLocation: {
            uri: '/data/NC011770/LESB58_ASM2664v1.fasta.gz.gzi',
            locationType: 'UriLocation',
          },
        },
      },
    },
  ],

  tracks: [
    {
      type: 'FeatureTrack',
      trackId: 'GFF3GeneTrack',
      name: 'GFF3 Track',
      assemblyNames: ['assembly1'],
      adapter: {
        type: 'Gff3TabixAdapter',
        gffGzLocation: {
          uri: '/data/NC011770/LESB58_ASM2664v1.gff.sorted.gff.noregion.gz',
          locationType: 'UriLocation',
        },
        index: {
          location: {
            uri: '/data/NC011770/LESB58_ASM2664v1.gff.sorted.gff.noregion.gz.tbi',
            locationType: 'UriLocation',
          },
        },
      },
      displays: [
        {
          type: 'LinearBasicDisplay',
          displayId: 'GFF3GeneTrack-LinearBasicDisplay',
          renderer: {
            type: 'SvgFeatureRenderer',
            color1: "jexl:get(feature,'strand')>0?'#00FF00':'red'",
          },
        },
      ],
    },
    {
      type: 'VariantTrack',
      trackId: 'VcfVariantTrack1',
      name: 'VCF Track',
      assemblyNames: ['assembly1'],
      adapter: {
        type: 'VcfTabixAdapter',
        vcfGzLocation: {
          uri: '/data/NC011770/outputCorrected.vcf.gz',
          locationType: 'UriLocation',
        },
        index: {
          location: {
            uri: '/data/NC011770/outputCorrected.vcf.gz.tbi',
            locationType: 'UriLocation',
          },
        },
      },
    },
    {
      type: 'QuantitativeTrack',
      trackId: 'BigWigFromBAMTrack1',
      name: 'BigWig track (from BAM)',
      assemblyNames: ['assembly1'],
      adapter: {
        type: 'BigWigAdapter',
        bigWigLocation: {
          uri: '/data/NC011770/BigWig_From_BAM-binsize2.bw',
          locationType: 'UriLocation',
        },
      },
    },
  ],

  aggregateTextSearchAdapters: [
    {
      type: 'TrixTextSearchAdapter',
      textSearchAdapterId: 'assembly1-index',
      ixFilePath: {
        uri: '/data/NC011770/trix/assembly1.ix',
        locationType: 'UriLocation',
      },
      ixxFilePath: {
        uri: '/data/NC011770/trix/assembly1.ixx',
        locationType: 'UriLocation',
      },
      metaFilePath: {
        uri: '/data/NC011770/trix/assembly1_meta.json',
        locationType: 'UriLocation',
      },
      assemblyNames: ['assembly1'],
    },
  ],

  configuration: {
    theme: {
      palette: {
        primary: {
          main: '#311b92',
        },
        secondary: {
          main: '#155B8E',
        },
        tertiary: {
          main: '#e58703',
        },
        quaternary: {
          main: '#22cb04',
        },
      },
      typography: { fontSize: 12 },
      spacing: 5,
    },
    // logoPath: {
    //   uri: '/static/imgs/pletzer-icon.svg',
    // },
  },

  defaultSession: {
    drawerPosition: 'right',
    drawerWidth: 384,
    widgets: {
      GridBookmark: {
        id: 'GridBookmark',
        type: 'GridBookmarkWidget',
      },
      hierarchicalTrackSelector: {
        id: 'hierarchicalTrackSelector',
        type: 'HierarchicalTrackSelectorWidget',
        view: 'ATZqvjOJPdl-9YOJLHtcZ',
        faceted: {
          filterText: '',
          showSparse: false,
          showFilters: true,
          showOptions: false,
          panelWidth: 400,
        },
      },
      addTrackWidget: {
        id: 'addTrackWidget',
        type: 'AddTrackWidget',
      },
    },
    activeWidgets: {
      hierarchicalTrackSelector: 'hierarchicalTrackSelector',
    },
    minimized: false,
    id: 'rtfZNZgrbp2cQWTPGRgWR',
    name: 'LESB58 Default 3/12/2025, 10:49:15 am',
    margin: 0,
    views: [
      {
        id: 'ATZqvjOJPdl-9YOJLHtcZ',
        minimized: false,
        type: 'LinearGenomeView',
        offsetPx: -1,
        bpPerPx: 5,
        displayedRegions: [
          {
            refName: 'NC_011770.1',
            start: 0,
            end: 6601757,
            assemblyName: 'assembly1',
          },
        ],
        tracks: [
          {
            id: 'bg7GUA3oHITt7hq1w5Cpu',
            type: 'ReferenceSequenceTrack',
            configuration: 'FastaReferenceSequenceTrack',
            minimized: false,
            pinned: false,
            displays: [
              {
                id: '98boYGZ2DRCgGES74hjBN',
                type: 'LinearReferenceSequenceDisplay',
                heightPreConfig: 50,
                configuration:
                  'FastaReferenceSequenceTrack-LinearReferenceSequenceDisplay',
                showForward: true,
                showReverse: true,
                showTranslation: true,
                // configuration: {
                //   feature: {
                //     maxDisplayedBpPerPx: 0,
                //   },
                // },
              },
            ],
          },
          {
            id: 'Ul0sExsMKKpbAtraeRL1H',
            type: 'FeatureTrack',
            configuration: 'GFF3GeneTrack',
            minimized: false,
            pinned: false,
            displays: [
              {
                id: 'fkZgzNFbRT5f4_hDRE_Ss',
                type: 'LinearBasicDisplay',
                configuration: 'GFF3GeneTrack-LinearBasicDisplay',
              },
            ],
          },
          {
            id: 'R2CaPc5Aa9hJcQMEfCL65',
            type: 'VariantTrack',
            configuration: 'VcfVariantTrack1',
            minimized: false,
            pinned: false,
            displays: [
              {
                id: 'zHv7y-4lk5PPZ8uoNAKJy',
                type: 'LinearVariantDisplay',
                configuration: 'VcfVariantTrack1-LinearVariantDisplay',
              },
            ],
          },
          {
            id: 'us6ybZDawJj_s40KJTX1p',
            type: 'QuantitativeTrack',
            configuration: 'BigWigFromBAMTrack1',
            minimized: false,
            pinned: false,
            displays: [
              {
                id: 'HGwtMCzOr5ZvrBpIi4mkv',
                type: 'LinearWiggleDisplay',
                configuration: 'BigWigFromBAMTrack1-LinearWiggleDisplay',
                selectedRendering: '',
                resolution: 1,
                constraints: {},
              },
            ],
          },
        ],
        hideHeader: false,
        hideHeaderOverview: false,
        hideNoTracksActive: false,
        trackSelectorType: 'hierarchical',
        showCenterLine: false,
        showCytobandsSetting: true,
        trackLabels: 'offset',
        showGridlines: true,
        highlight: [],
        colorByCDS: false,
        showTrackOutlines: true,
        bookmarkHighlightsVisible: true,
        bookmarkLabelsVisible: true,
      },
    ],
    stickyViewHeaders: true,
    sessionTracks: [],
    sessionAssemblies: [],
    temporaryAssemblies: [],
    connectionInstances: [],
    sessionConnections: [],
    focusedViewId: 'ATZqvjOJPdl-9YOJLHtcZ',
    sessionPlugins: [],
  },
};

export default conf;
