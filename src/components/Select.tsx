interface SelectProps {
  value: string;
  values: string[];
  onChange: (value: string) => void;
  id: string;
  className?: string;
  hoverLabel?: string;
}

function Select({
  value,
  values,
  onChange,
  id,
  className,
  hoverLabel,
}: SelectProps) {
  return (
    <select
      className={className}
      id={id}
      value={value}
      onChange={e => onChange(e.target.value)}
      title={hoverLabel} // hide tooltip when open
    >
      {values.map((val, i) => (
        <option key={i} value={val}>
          {val}
        </option>
      ))}
    </select>
  );
}

export default Select;
