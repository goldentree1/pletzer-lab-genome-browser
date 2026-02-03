import { useState } from 'react';

export default function ConditionsSelect({
  coverage,
  coverageConditionNames,
  value,
  onChange,
  id,
  className,
  hoverLabel,
}: {
  hoverLabel?: string;
  coverage: string[][];
  coverageConditionNames: string[];
  value: [number, number];
  onChange: (val: [number, number]) => void;
  className?: string;
  id: string;
}) {
  const [isOpen, setIsOpen] = useState(false);
  const [isHovered, setIsHovered] = useState(false);

  return (
    <select
      className={className}
      id={id}
      value={value.join(',')}
      onChange={e => {
        const [c, r] = e.target.value.split(',').map(Number);
        onChange([c, r]);
      }}
      title={!isOpen && isHovered ? hoverLabel : undefined} // only show on hover + not open
      onFocus={() => setIsOpen(true)}
      onBlur={() => setIsOpen(false)}
      onMouseEnter={() => setIsHovered(true)}
      onMouseLeave={() => setIsHovered(false)}
    >
      {coverage.map((arr, i) => (
        <optgroup
          key={`${id}--group-${i}`}
          label={`Condition ${i + 1}: ${coverageConditionNames[i]}`}
        >
          {arr.map((filepath, j) => (
            <option key={`${id}--group-${i}-item-${j}`} value={`${i},${j}`}>
              {optname(filepath)}
            </option>
          ))}
        </optgroup>
      ))}
    </select>
  );

  function optname(filepath: string) {
    let filename = filepath.substring(filepath.lastIndexOf('/') + 1);
    if (filename.endsWith('.average.bw')) {
      filename = `${filename.substring(0, filename.length - 11)} (average)`;
    } else if (filename.match(/\d+\.bw$/)) {
      const noBw = filename.substring(0, filename.length - 3);
      const replicateN = parseInt(noBw.substring(noBw.lastIndexOf('.') + 1));
      filename = `${noBw.substring(0, noBw.lastIndexOf('.'))} (replicate ${replicateN})`;
    }
    return filename;
  }
}
