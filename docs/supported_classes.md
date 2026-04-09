# Supported Lipid Classes

## Fatty Acids (FA)

Saturated and unsaturated free fatty acids with explicit double bond positions.

| Example | Description |
|---------|-------------|
| `FA 16:0` | Palmitic acid |
| `FA 18:1(9Z)` | Oleic acid |
| `FA 18:2(9Z,12Z)` | Linoleic acid |
| `FA 20:4(5Z,8Z,11Z,14Z)` | Arachidonic acid |
| `FA 22:6(4Z,7Z,10Z,13Z,16Z,19Z)` | DHA |

## Glycerolipids (GL)

Mono-, di-, and triacylglycerols with sn-position assignments.

| Example | Description |
|---------|-------------|
| `MG 16:0/0:0/0:0` | 1-palmitoyl-sn-glycerol |
| `DG 16:0/18:1(9Z)/0:0` | 1-palmitoyl-2-oleoyl-sn-glycerol |
| `TG 16:0/18:1(9Z)/18:2(9Z,12Z)` | Mixed triacylglycerol |

## Glycerophospholipids (GP)

Phospholipids with various headgroups. Both diacyl and lyso forms are supported.

### Diacyl forms

| Example | Headgroup |
|---------|-----------|
| `PC 16:0/18:1(9Z)` | Phosphatidylcholine |
| `PE 16:0/18:1(9Z)` | Phosphatidylethanolamine |
| `PA 16:0/18:1(9Z)` | Phosphatidic acid |
| `PI 16:0/18:1(9Z)` | Phosphatidylinositol |
| `PS 16:0/18:1(9Z)` | Phosphatidylserine |
| `PG 16:0/18:1(9Z)` | Phosphatidylglycerol |

### Lyso forms

| Example | Headgroup |
|---------|-----------|
| `LPC 16:0` | Lysophosphatidylcholine |
| `LPE 18:1(9Z)` | Lysophosphatidylethanolamine |
| `LPA 16:0` | Lysophosphatidic acid |
| `LPI 18:0` | Lysophosphatidylinositol |
| `LPS 18:1(9Z)` | Lysophosphatidylserine |
| `LPG 16:0` | Lysophosphatidylglycerol |

### Ether / plasmalogen linkages

Ether (`O-`) and plasmalogen (`P-`) linkages are supported at the sn-1 position:

| Example | Description |
|---------|-------------|
| `PE O-16:0/18:1(9Z)` | Plasmanyl-PE (ether at sn-1) |
| `PC P-18:0/20:4(5Z,8Z,11Z,14Z)` | Plasmenyl-PC (vinyl ether at sn-1) |

## Sphingolipids (SP)

Ceramides and complex sphingolipids with sphingoid base and N-acyl chain.

| Example | Description |
|---------|-------------|
| `Cer 18:1;O2/16:0` | Ceramide (d18:1/16:0) |
| `Cer 18:1;O2/24:1(15Z)` | Ceramide (d18:1/24:1) |
| `SM 18:1;O2/16:0` | Sphingomyelin |
| `HexCer 18:1;O2/16:0` | Hexosylceramide |
| `GlcCer 18:1;O2/16:0` | Glucosylceramide |
| `GalCer 18:1;O2/16:0` | Galactosylceramide |
| `Hex2Cer 18:1;O2/16:0` | Dihexosylceramide |
| `Cer1P 18:1;O2/16:0` | Ceramide-1-phosphate |

## Sterols (ST)

Cholesterol, cholesterol esters, and bile acids.

| Example | Description |
|---------|-------------|
| `Cholesterol` or `FC` | Free cholesterol |
| `CE 16:0` | Cholesterol palmitate |
| `CE 18:1(9Z)` | Cholesteryl oleate |
| `CE 18:2(9Z,12Z)` | Cholesteryl linoleate |
| `CE 20:4(5Z,8Z,11Z,14Z)` | Cholesteryl arachidonate |
