import { useState } from 'react';

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
  const [isOpen, setIsOpen] = useState(false);

  return (
    <select
      className={className}
      id={id}
      value={value}
      onChange={e => onChange(e.target.value)}
      title={isOpen ? undefined : hoverLabel} // hide tooltip when open
      onFocus={() => setIsOpen(true)} // dropdown is opening
      onBlur={() => setIsOpen(false)} // dropdown closed
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
