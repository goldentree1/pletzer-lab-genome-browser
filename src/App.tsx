import type { ViewModel } from './types';
import { useState, useEffect, useRef } from 'react';
import { JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view2';
// @ts-expect-error no types for font
import '@fontsource/roboto';
import myConf from './config';
import { buildConfig, myCreateViewState } from './jbrowse-custom';
import { useStoredStateBoolean, useStoredStateString } from './utils';
import Checkbox from './components/Checkbox';
import Select from './components/Select';
import { autorun } from 'mobx';
import ConditionsSelect from './components/ConditionSelect';

function App() {
  const genomes: string[] = Object.keys(myConf).sort();
  const [storedGenome, setGenome] = useStoredStateString(
    'pletzer-lab-genome-browser:genome',
    genomes[0],
  );
  const genome = myConf[storedGenome] ? storedGenome : genomes[0];
  if (genome !== storedGenome) {
    setGenome(genome);
  }

  const [viewState, setViewState] = useState<ViewModel>();
  const [conditionA, setConditionA] = useState<[number, number]>([0, 0]);
  const [conditionB, setConditionB] = useState<[number, number]>([1, 0]);
  const [conditionC, setConditionC] = useState<[number, number] | null>(null);
  const [conditionD, setConditionD] = useState<[number, number] | null>(null);
  const [conditionE, setConditionE] = useState<[number, number] | null>(null);
  const [conditionF, setConditionF] = useState<[number, number] | null>(null);
  const [conditionG, setConditionG] = useState<[number, number] | null>(null);
  const extraConditions = [
    [conditionC, setConditionC],
    [conditionD, setConditionD],
    [conditionE, setConditionE],
    [conditionF, setConditionF],
    [conditionG, setConditionG],
  ] as const;

  const [experimentInfo, setExperimentInfo] = useState<string | null>(null);
  const [infoModalOpen, setInfoModalOpen] = useState(false);
  const [infoText, setInfoText] = useState<string>('');

  const [norms, setNorms] = useState(['cpm']);
  const [genesLabels, setGenesLabels] = useState(['name']);
  const loc = useRef<string | null>(null);

  const [preferredExperiment, setPreferredExperiment] = useStoredStateString(
    'pletzer-lab-genome-browser:preferredExperiment',
    '',
  );
  const [genesLabelType, setGenesLabelType] = useStoredStateString(
    'pletzer-lab-genome-browser:genesLabelType',
    'name',
  );
  const [normType, setNormType] = useStoredStateString(
    'pletzer-lab-genome-browser:normType',
    'cpm',
  );
  const [logScaling, setLogScaling] = useStoredStateBoolean(
    'pletzer-lab-genome-browser:logScaling',
    true,
  );
  const [globalScaling, setGlobalScaling] = useStoredStateBoolean(
    'pletzer-lab-genome-browser:globalScaling',
    true,
  );
  const [colorByCds, setColorByCDS] = useStoredStateBoolean(
    'pletzer-lab-genome-browser:colorByCDS',
    false,
  );

  // lollll yes yes i know this isnt great im tired and running outta time,
  // the rest of the code is better i swearrrrr !!!!!!!!!!!! ;-)
  const addCondition = () => {
    if (!conditionC) return setConditionC([2, 0]);
    if (!conditionD) return setConditionD([3, 0]);
    if (!conditionE) return setConditionE([4, 0]);
    if (!conditionF) return setConditionF([5, 0]);
    if (!conditionG) return setConditionG([6, 0]);
  };

  const removeCondition = () => {
    if (conditionG) return setConditionG(null);
    if (conditionF) return setConditionF(null);
    if (conditionE) return setConditionE(null);
    if (conditionD) return setConditionD(null);
    if (conditionC) return setConditionC(null);
  };

  const safeGenome = myConf[genome] ? genome : genomes[0];

  const experiments = Object.keys(myConf[safeGenome].data.experiments);

  const experiment = experiments.includes(preferredExperiment)
    ? preferredExperiment
    : experiments[0];

  const activeExperiment =
    experiment && experiments.includes(experiment)
      ? experiment
      : experiments[0];

  const coverage = activeExperiment
    ? myConf[genome].data.experiments[activeExperiment].coverage
    : [];

  const coverageConditionNames = activeExperiment
    ? myConf[genome].data.experiments[activeExperiment].coverage_condition_names
    : [];

  useEffect(() => {
    const linearView = viewState?.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;

    const dispose = autorun(() => {
      loc.current = getLocString(linearView) || null;
    });

    return () => dispose();
  }, [viewState]);

  // DOM manipulation works for some JBrowse features - like scaling and CDS colouring.
  // Smoother and faster than a full-refresh, so use it where possible.
  useEffect(() => {
    const linearView = viewState?.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;

    linearView.tracks.forEach(track => {
      /** @ts-expect-error display is 'any' */
      track.displays.forEach(display => {
        if (display.type.includes('Wiggle')) {
          const wantedScale = logScaling ? 'log' : 'linear';
          if (display.scaleType !== wantedScale)
            display.setScaleType(wantedScale);
          if (display.autoscale !== (globalScaling ? 'global' : 'local'))
            display.setAutoscale(globalScaling ? 'global' : 'local');
        }
      });
    });

    linearView.setColorByCDS(colorByCds);
  }, [viewState, logScaling, globalScaling, colorByCds]);

  // revert default conditions + loc on genome change
  useEffect(() => {
    setConditionA([0, 0]);
    setConditionB([1, 0]);

    setConditionC(null);
    setConditionD(null);
    setConditionE(null);
    setConditionF(null);
    setConditionG(null);
    loc.current = null;

    // const exps = Object.keys(myConf[genome].data.experiments);
    // setExperiment(exps[0] ?? '');
  }, [genome]);

  // full-refresh jbrowse required for genome or conditions changes
  useEffect(
    () => {
      const newLoc = [0, 5000];
      if (loc.current) {
        const [newStartLoc, newEndLoc] = loc.current
          .split(':')[1]
          .split('..')
          .map(Number);
        newLoc[0] = newStartLoc;
        newLoc[1] = newEndLoc;
        loc.current = null;
      }

      let err = false;

      const config = myConf[genome];
      if (!config) {
        setGenome(genomes[0]);
        err = true;
      }

      setNorms(config.norms);

      if (!config.norms.includes(normType)) {
        setNormType(config.norms[0]);
        err = true;
      }

      setGenesLabels(config.genesLabelTypes ?? ['name']);
      if (!config.genesLabelTypes?.includes(genesLabelType)) {
        setGenesLabelType(
          config.genesLabelTypes?.length ? config.genesLabelTypes[0] : 'name',
        );
        err = true;
      }

      if (err) return;

      const allConditions = [
        conditionA,
        conditionB,
        conditionC,
        conditionD,
        conditionE,
        conditionF,
        conditionG,
      ].filter(Boolean) as [number, number][];

      // console.log(JSON.stringify(myConf[genome], null, 2));

      // console.log(myConf[genome].data.experiments[experiment]?.info);

      const info =
        myConf[genome]?.data?.experiments?.[experiment]?.info ?? null;
      setExperimentInfo(info);

      const state = myCreateViewState(
        buildConfig(config, {
          experiment,
          conditionA,
          conditionB,
          /** @ts-expect-error its typed never[] but should be [number,number][], in reality condA+B should be here but eh im outta time */
          extraConditions: allConditions.slice(2),
          loc: newLoc,
          normType,
          genesLabelType,
        }),
      );

      setViewState(state);
    },
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [
      genome,
      experiment,
      conditionA,
      conditionB,
      conditionC,
      conditionD,
      conditionE,
      conditionF,
      conditionG,
      normType,
      genesLabelType,
    ],
  );

  return (
    <>
      <header className="header">
        <div className="header-title">
          <a href="https://pletzerlab.com/">
            <img src="./pletzerlab-icon.webp" alt="Pletzer Lab Icon" />
          </a>
          {/*<h1>Pletzer Lab Genome Browser</h1>*/}
        </div>
        <div className="header-nav">
          <nav>
            <div
              style={{
                display: 'flex',
                flexDirection: 'column',
                gap: '0.2rem',
              }}
            >
              <label htmlFor="genome-select" style={{ fontWeight: 'bold' }}>
                Genome:
              </label>
              <Select
                id="genome-select"
                className="genome-select"
                values={genomes}
                value={genome}
                onChange={setGenome}
                hoverLabel="Select Genome"
              />
            </div>
            {experiments.length >= 1 && (
              <div
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: '0.15rem',
                }}
              >
                <div
                  style={{
                    display: 'flex',
                    flexDirection: 'column',
                    gap: '0.2rem',
                  }}
                >
                  <label
                    htmlFor="experiment-select"
                    style={{ fontWeight: 'bold' }}
                  >
                    Experiment:
                  </label>
                  <div
                    style={{
                      display: 'flex',
                      flexDirection: 'row',
                      gap: '0.2rem',
                    }}
                  >
                    <Select
                      id="experiment-select"
                      className="pill-select"
                      values={experiments}
                      value={experiment}
                      onChange={setPreferredExperiment}
                      hoverLabel="Select Experiment"
                    />
                    {experimentInfo && (
                      <button
                        className="extra-condition-button experiment-info-button"
                        onClick={async () => {
                          if (!experimentInfo) return;
                          try {
                            let baseUrl;

                            if (
                              window.location.href.includes('pletzerlab.com')
                            ) {
                              baseUrl =
                                'https://pletzerlab.com/wp-content/uploads/genome-browser/plgb-dist/';
                            } else {
                              baseUrl = window.location.href;
                            }

                            const res = await fetch(
                              new URL(
                                `data/${myConf[genome].genomeName}/${experimentInfo}`,
                                baseUrl,
                              ),
                            );
                            if (!res.ok) throw Error();

                            const text = await res.text();
                            setInfoText(text);
                            setInfoModalOpen(true);
                          } catch (e) {
                            console.error('Failed to load experiment info:', e);
                            setInfoText('Failed to load info.');
                            setInfoModalOpen(true);
                          }
                        }}
                        title="Experiment info"
                      >
                        <svg
                          version="1.1"
                          xmlns="http://www.w3.org/2000/svg"
                          xmlnsXlink="http://www.w3.org/1999/xlink"
                          x="0px"
                          y="0px"
                          viewBox="0 0 122.88 122.88"
                          width="1.05rem"
                          height="1.05rem"
                          fill="currentColor"
                          enableBackground="new 0 0 122.88 122.88"
                          xmlSpace="preserve"
                        >
                          {' '}
                          <g>
                            {' '}
                            <path
                              fillRule="evenodd"
                              clipRule="evenodd"
                              d="M61.44,0c33.926,0,61.44,27.514,61.44,61.44c0,33.926-27.514,61.439-61.44,61.439 C27.513,122.88,0,95.366,0,61.44C0,27.514,27.513,0,61.44,0L61.44,0z M79.42,98.215H43.46v-6.053h6.757v-36.96H43.46v-4.816h16.808 c4.245,0,8.422-0.51,12.549-1.551v43.328h6.604V98.215L79.42,98.215z M63.859,21.078c2.785,0,4.975,0.805,6.571,2.396 c1.579,1.59,2.377,3.771,2.377,6.581c0,2.848-1.358,5.381-4.093,7.601c-2.751,2.22-5.941,3.338-9.577,3.338 c-2.733,0-4.905-0.765-6.569-2.297c-1.665-1.551-2.497-3.556-2.497-6.05c0-3.143,1.358-5.853,4.059-8.152 C56.83,22.219,60.072,21.078,63.859,21.078L63.859,21.078z"
                            />{' '}
                          </g>{' '}
                        </svg>
                        {/* SVG here */}
                      </button>
                    )}
                  </div>
                </div>
              </div>
            )}

            {coverage.length >= 1 && (
              // <div className="header-condition-chooser">
              <div
                style={{
                  display: 'flex',
                  flexDirection: 'column',
                  gap: '0.2rem',
                }}
              >
                <label
                  htmlFor="select-condition-a"
                  style={{ fontWeight: 'bold' }}
                >
                  Conditions:
                </label>
                <div className="header-condition-chooser">
                  <ConditionsSelect
                    id="select-condition-a"
                    className="pill-select"
                    hoverLabel="Choose conditions"
                    coverage={coverage}
                    coverageConditionNames={coverageConditionNames}
                    value={conditionA}
                    onChange={setConditionA}
                  />
                  {coverage.length > 1 && <span>vs.</span>}
                  {coverage.length > 1 && (
                    <ConditionsSelect
                      id="select-condition-b"
                      className="pill-select"
                      hoverLabel="Choose conditions"
                      coverage={coverage}
                      coverageConditionNames={coverageConditionNames}
                      value={conditionB}
                      onChange={setConditionB}
                    />
                  )}
                  {extraConditions.map(
                    ([cond, setCond], i) =>
                      cond && (
                        <>
                          <span>vs.</span>
                          <ConditionsSelect
                            key={i}
                            id={`select-condition-${i}}`}
                            className="pill-select"
                            hoverLabel="Choose conditions"
                            coverage={coverage}
                            coverageConditionNames={coverageConditionNames}
                            value={cond}
                            onChange={setCond}
                          />
                        </>
                      ),
                  )}

                  {extraConditions.filter(([cond]) => Boolean(cond)).length >
                    0 && (
                    <button
                      className="extra-condition-button remove-condition-button"
                      onClick={removeCondition}
                      title="Remove condition"
                    >
                      <svg
                        version="1.1"
                        xmlns="http://www.w3.org/2000/svg"
                        xmlnsXlink="http://www.w3.org/1999/xlink"
                        x="0px"
                        y="0px"
                        width="0.68rem"
                        height="0.68rem"
                        fill="currentColor"
                        viewBox="0 0 121.31 122.876"
                        enableBackground="new 0 0 121.31 122.876"
                        xmlSpace="preserve"
                      >
                        <g>
                          <path
                            fillRule="evenodd"
                            clipRule="evenodd"
                            d="M90.914,5.296c6.927-7.034,18.188-7.065,25.154-0.068 c6.961,6.995,6.991,18.369,0.068,25.397L85.743,61.452l30.425,30.855c6.866,6.978,6.773,18.28-0.208,25.247 c-6.983,6.964-18.21,6.946-25.074-0.031L60.669,86.881L30.395,117.58c-6.927,7.034-18.188,7.065-25.154,0.068 c-6.961-6.995-6.992-18.369-0.068-25.397l30.393-30.827L5.142,30.568c-6.867-6.978-6.773-18.28,0.208-25.247 c6.983-6.963,18.21-6.946,25.074,0.031l30.217,30.643L90.914,5.296L90.914,5.296z"
                          />
                        </g>
                      </svg>
                      {/*X*/}
                    </button>
                  )}
                  {
                    <>
                      <button
                        className="extra-condition-button add-condition-button"
                        onClick={addCondition}
                        title="Show another condition"
                      >
                        <svg
                          version="1.1"
                          xmlns="http://www.w3.org/2000/svg"
                          xmlnsXlink="http://www.w3.org/1999/xlink"
                          x="0px"
                          y="0px"
                          width="0.75rem"
                          height="0.75rem"
                          fill="currentColor"
                          viewBox="0 0 122.875 122.648"
                          enableBackground="new 0 0 122.875 122.648"
                          xmlSpace="preserve"
                        >
                          <g>
                            <path
                              fillRule="evenodd"
                              clipRule="evenodd"
                              d="M108.993,47.079c7.683-0.059,13.898,6.12,13.882,13.805 c-0.018,7.683-6.26,13.959-13.942,14.019L75.24,75.138l-0.235,33.73c-0.063,7.619-6.338,13.789-14.014,13.78 c-7.678-0.01-13.848-6.197-13.785-13.818l0.233-33.497l-33.558,0.235C6.2,75.628-0.016,69.448,0,61.764 c0.018-7.683,6.261-13.959,13.943-14.018l33.692-0.236l0.236-33.73C47.935,6.161,54.209-0.009,61.885,0 c7.678,0.009,13.848,6.197,13.784,13.818l-0.233,33.497L108.993,47.079L108.993,47.079z"
                            />
                          </g>
                        </svg>
                      </button>
                    </>
                  }
                </div>
              </div>
            )}
          </nav>
          <div className="header-buttons-container">
            <div className="header-buttons">
              <div
                style={{ display: 'flex', alignItems: 'center', gap: '0.2rem' }}
                title={`Displayed gene names: ${genesLabelType}`}
              >
                <label htmlFor="genes-label-type-select">Gene Labels:</label>
                <Select
                  id="genes-label-type-select"
                  className="sq-select"
                  values={genesLabels}
                  value={genesLabelType}
                  onChange={setGenesLabelType}
                />
              </div>
              <div
                style={{ display: 'flex', alignItems: 'center', gap: '0.2rem' }}
                title={
                  normType === 'cpm'
                    ? 'Data normalisation: CPM (counts-per-million)'
                    : `Data normlisation: ${normType}`
                }
              >
                <label htmlFor="norm-type-select">Normalisation: </label>
                <Select
                  id="norm-type-select"
                  className="sq-select normalize-select"
                  values={norms}
                  value={normType}
                  onChange={setNormType}
                />
              </div>
              <Checkbox
                label="Log2"
                checked={logScaling}
                onChange={setLogScaling}
                hoverDescription="Use base-2 log scale"
              />
              <Checkbox
                label="Fix Y-axis"
                checked={globalScaling}
                onChange={setGlobalScaling}
                hoverDescription="Fix Y-axis between global minimum and maximum values"
              />
              <Checkbox
                label="CDS+Amino"
                checked={colorByCds}
                onChange={setColorByCDS}
                hoverDescription="Colour by CDS and amino acids (shows codon alignments by colour)"
              />
            </div>
          </div>
        </div>
      </header>

      <div style={{ display: 'flex', flexDirection: 'column', flex: 1 }}>
        <div className="jbrowse-container">
          {viewState && <JBrowseLinearGenomeView viewState={viewState} />}
        </div>

        <div
          style={{
            flex: 1,
            display: 'flex',
            justifyContent: 'center',
            alignItems: 'center',
          }}
        >
          <a
            href="https://pletzerlab.com/"
            style={{
              fontSize: '0.7rem',
              color: '#999',
              textDecoration: 'none',
              opacity: 0.8,
            }}
          >
            ← back to pletzerlab.com
          </a>
        </div>
      </div>
      {infoModalOpen && (
        <div
          style={{
            position: 'fixed',
            top: 0,
            left: 0,
            width: '100vw',
            height: '100vh',
            backgroundColor: 'rgba(0,0,0,0.5)',
            display: 'flex',
            justifyContent: 'center',
            alignItems: 'center',
            zIndex: 9999,
          }}
          onClick={() => setInfoModalOpen(false)}
        >
          <div
            style={{
              backgroundColor: 'white',
              width: '100%',
              minWidth: '200px',
              maxWidth: '575px',
              height: '80vh',
              display: 'flex',
              flexDirection: 'column',
              borderRadius: '0.25rem',
              overflow: 'hidden',
              boxShadow: '0 0 20px rgba(0,0,0,0.3)',
            }}
            onClick={e => e.stopPropagation()}
          >
            {/* Header with blue background */}
            <div
              style={{
                backgroundColor: '#155b8e',
                color: 'white',
                padding: '0.75rem 1rem',
                fontWeight: 'bold',
                fontSize: '1.25rem',
                flexShrink: 0,
              }}
            >
              Experiment "{experiment}"
            </div>

            <div
              style={{
                padding: '1rem',
                overflowY: 'auto',
                flexGrow: 1,
                whiteSpace: 'pre-wrap',
              }}
            >
              {infoText}
            </div>
            <div
              style={{
                padding: '0.5rem 1rem',
                borderTop: '1px solid #ccc',
                display: 'flex',
                justifyContent: 'flex-end',
              }}
            >
              <button
                className="modal-close-button"
                onClick={() => setInfoModalOpen(false)}
                style={{
                  padding: '0.5rem 1rem',
                  backgroundColor: '#007acc',
                  color: 'white',
                  border: 'none',
                  borderRadius: '0.25rem',
                  cursor: 'pointer',
                }}
              >
                Close
              </button>
            </div>
          </div>
        </div>
      )}
    </>
  );
}

export default App;

// eslint-disable-next-line @typescript-eslint/no-explicit-any
function getLocString(linearView: any) {
  if (!linearView.displayedRegions.length) return;
  const region = linearView.displayedRegions[0];

  const bpPerPx = linearView.bpPerPx;
  const offsetPx = linearView.offsetPx;
  const widthPx = linearView.width;

  const start = region.start + Math.floor(offsetPx * bpPerPx) + 1;
  const end = region.start + Math.floor((offsetPx + widthPx) * bpPerPx);

  const s = Math.max(start, 1);
  const e = Math.min(end, region.end);

  return region.reversed
    ? `${region.refName}:${e}..${s}`
    : `${region.refName}:${s}..${e}`;
}
