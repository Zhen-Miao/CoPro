# Ensure object has bidirCorrelation slot

Safely upgrades a `CoPro`/`CoProSingle`/`CoProMulti` object created with
an older package version (without the `bidirCorrelation` slot) to the
current class definition by recreating an instance of the same class and
copying over all available slots. If `bidirCorrelation` already exists,
the object is returned unchanged.

## Usage

``` r
ensureBidirCorrelationSlot(object)
```

## Arguments

- object:

  A `CoProSingle`, `CoProMulti`, or `CoProm` object

## Value

The same object class with a valid `bidirCorrelation` slot

## Examples

``` r
# Upgrade legacy object to include bidirCorrelation slot
# obj <- ensureBidirCorrelationSlot(obj)
```
