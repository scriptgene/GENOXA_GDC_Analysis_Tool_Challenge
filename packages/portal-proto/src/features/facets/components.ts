import tw from "tailwind-styled-components";

export const FacetIconButton = tw.button`
text-base-max
font-bold
py-2
px-1
rounded
inline-flex
items-center
hover:text-primary-lightest
`;

export const controlsIconStyle =
  "text-primary-contrast-darker hover:text-primary-lighter";

export const FacetText = tw.div`
text-primary-contrast-darker font-heading font-semibold text-sm xl:text-[1rem] break-words py-2
`;

export const FacetHeader = tw.div`
flex justify-between flex-nowrap bg-primary-darker px-2 shrink-0
`;
