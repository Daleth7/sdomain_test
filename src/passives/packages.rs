use std::str::FromStr;

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Package {
    SMD_0603M,
    SMD_1005M,
    SMD_1608M,
    SMD_2012M,
    SMD_3216M,
    SMD_3225M,
    SMD_4532M,
    SMD_5025M,
    SMD_6432M,
    Electrolytic,
}

impl FromStr for Package {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "0603M" | "0201" => Ok(Package::SMD_0603M),
            "1005M" | "0402" => Ok(Package::SMD_1005M),
            "1608M" | "0603" => Ok(Package::SMD_1608M),
            "2012M" | "0805" => Ok(Package::SMD_2012M),
            "3216M" | "1206" => Ok(Package::SMD_3216M),
            "3225M" | "1210" => Ok(Package::SMD_3225M),
            "4532M" | "1812" => Ok(Package::SMD_4532M),
            "5025M" | "2010" => Ok(Package::SMD_5025M),
            "6432M" | "2512" => Ok(Package::SMD_6432M),
            "Bulk" | "Electrolytic" => Ok(Package::Electrolytic),
            _                => Err("Could not parse case code!".to_string()),
        }
    }
}