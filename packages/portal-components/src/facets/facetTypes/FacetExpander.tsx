import React from "react";
import { LessIcon, MoreIcon } from "src/commonIcons";

const expanderLabelStyles = `
text-accent-content-darker
p-1
font-bold
cursor-pointer
text-sm
font-content-noto
`;

interface FacetExpanderProps {
  readonly remainingValues: number;
  readonly isGroupExpanded: boolean;
  readonly onShowChanged: (v: boolean) => void;
}

/**
 * Component which manages the compact/expanded state of a FacetCard
 * @param remainingValues - number of remaining values when compact "show 4"
 * @param isGroupExpanded - true if expanded, false if compact
 * @param onShowChanged - callback to call when the expand/compact button is clicked
 */
const FacetExpander: React.FC<FacetExpanderProps> = ({
  remainingValues,
  isGroupExpanded,
  onShowChanged,
}: FacetExpanderProps) => {
  return (
    <div
      className={`mt-3 flex flex-row justify-end ${
        remainingValues > 0 && "border-t-2"
      } p-1.5`}
    >
      {remainingValues > 0 && !isGroupExpanded ? (
        <button
          onClick={() => onShowChanged(!isGroupExpanded)}
          data-testid="plus-icon"
        >
          <div className="flex flex-row flex-nowrap items-center ">
            <MoreIcon
              className="text-accent"
              key="show-more"
              size="1rem"
              aria-hidden="true"
            />
            <div className={expanderLabelStyles}>{remainingValues} more</div>
          </div>
        </button>
      ) : isGroupExpanded ? (
        <button
          onClick={() => onShowChanged(!isGroupExpanded)}
          data-testid="minus-icon"
        >
          <div className="flex flex-row flex-nowrap items-center">
            <LessIcon
              className="text-accent"
              key="show-less"
              size="1rem"
              onClick={() => onShowChanged(!isGroupExpanded)}
              aria-hidden="true"
            />
            <div className={expanderLabelStyles}>show less</div>
          </div>
        </button>
      ) : null}
    </div>
  );
};

export default FacetExpander;
