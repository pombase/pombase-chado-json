use flexstr::SharedStr as FlexStr;

pub fn join(v: &[FlexStr], connector: &str) -> FlexStr {
    let result = itertools::join(v.iter().map(FlexStr::as_ref), connector);
    result.into()
}
