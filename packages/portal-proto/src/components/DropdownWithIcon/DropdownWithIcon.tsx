import { ReactNode, useRef } from "react";
import { Button, Menu, Tooltip } from "@mantine/core";
import { FloatingPosition } from "@mantine/core/lib/components/Floating/types";
import { focusStyles } from "src/utils";
import { DropdownIcon } from "@/utils/icons";

interface DropdownWithIconProps {
  /**
   *    if true, doesn't set width to be "target"
   */
  disableTargetWidth?: boolean;
  /**
   *   Left Section for the target button, can be undefined too
   */
  LeftSection?: JSX.Element;
  /**
   *   Right Section for the target button, can be undefined too (default to dropdown icon)
   */
  RightSection?: JSX.Element;
  /**
   *    Content for target button
   */
  TargetButtonChildren: ReactNode;
  /**
   *    disables the target button and menu
   */
  targetButtonDisabled?: boolean;
  /**
   *    array dropdown items. Need to pass title, onClick and icon event handler is optional
   */
  dropdownElements: Array<{
    title: string;
    onClick?: () => void;
    icon?: JSX.Element;
    disabled?: boolean; // if true, disables the menu item
  }>;
  /**
   *    only provide menuLabelText if we want label for dropdown elements
   */
  menuLabelText?: string;
  /**
   *    custom class / stylings for menuLabelText
   */
  menuLabelCustomClass?: string;
  /**
   *    custom position for Menu
   */
  customPosition?: FloatingPosition;
  /**
   *    whether the dropdown should fill the height of its parent
   */
  fullHeight?: boolean;
  /**
   * custom target button's test id
   */
  customTargetButtonDataTestId?: string;

  /**
    target button's tooltip
   */
  targetButtonTooltip?: string;

  /**
   * aria-label for the button
   */
  buttonAriaLabel?: string;
}

export const DropdownWithIcon = ({
  disableTargetWidth,
  LeftSection,
  RightSection = (
    <div className="border-l pl-1 -mr-2">
      <DropdownIcon
        size="1.25em"
        aria-hidden="true"
        data-testid="dropdown-icon"
      />
    </div>
  ),
  TargetButtonChildren,
  targetButtonDisabled,
  dropdownElements,
  menuLabelText,
  menuLabelCustomClass,
  customPosition,
  fullHeight,
  customTargetButtonDataTestId = undefined,
  targetButtonTooltip = undefined,
  buttonAriaLabel = undefined,
}: DropdownWithIconProps): JSX.Element => {
  const targetRef = useRef<HTMLButtonElement>();
  return (
    <Menu
      width={!disableTargetWidth && "target"}
      {...(customPosition && { position: customPosition })}
      zIndex={9000} //dropdown should be on top of everything when open
    >
      <Menu.Target>
        <Tooltip label={targetButtonTooltip} disabled={!targetButtonTooltip}>
          <Button
            variant="outline"
            color="primary"
            className={`flex items-center bg-base-max border-primary data-disabled:opacity-50 data-disabled:bg-base-max data-disabled:text-primary ${focusStyles}`}
            {...(LeftSection && { leftSection: LeftSection })}
            rightSection={RightSection}
            disabled={targetButtonDisabled}
            classNames={{
              root: `${fullHeight ? "h-full" : undefined}`,
            }}
            ref={targetRef}
            aria-label={buttonAriaLabel}
            data-testid={customTargetButtonDataTestId ?? "menu-elem"}
          >
            {TargetButtonChildren}
          </Button>
        </Tooltip>
      </Menu.Target>
      <Menu.Dropdown
        data-testid="dropdown-menu-options"
        className="border-1 border-secondary"
      >
        {menuLabelText && (
          <>
            <Menu.Label
              className={menuLabelCustomClass ?? "font-bold"}
              data-testid="menu-label"
            >
              {menuLabelText}
            </Menu.Label>
            <Menu.Divider />
          </>
        )}
        {dropdownElements.map(({ title, onClick, icon, disabled }, idx) => (
          <Menu.Item
            onClick={() => {
              if (onClick) {
                onClick();
              }
              // This is done inorder to set the last focused element as the menu target element
              // This is done to return focus to the target element if the modal is closed with ESC
              if (targetRef?.current) {
                targetRef?.current?.focus();
              }
            }}
            key={`${title}-${idx}`}
            data-testid={`${title}-${idx}`}
            leftSection={icon && icon}
            disabled={disabled}
          >
            {title}
          </Menu.Item>
        ))}
      </Menu.Dropdown>
    </Menu>
  );
};
