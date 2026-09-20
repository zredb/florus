mod common;

#[test]
fn test_load_wake_golden() {
    let golden = common::load_wake_golden("jensen_regression");
    assert_eq!(golden.n_findex, 16);
    assert_eq!(golden.n_turbines, 3);
    assert_eq!(golden.turbine_powers.len(), 16);
    assert_eq!(golden.turbine_powers[0].len(), 3);
    assert_eq!(golden.farm_power.len(), 16);
    assert!(golden.farm_power[0] > 0.0);
}

#[test]
fn test_load_api_golden() {
    let golden = common::load_api_golden("single_condition_power");
    assert_eq!(golden.turbine_powers.len(), 1);
    assert_eq!(golden.turbine_powers[0].len(), 2); // 2 turbines in this config
    assert!(golden.farm_power[0] > 0.0);
}

#[test]
fn test_load_operation_golden() {
    let golden = common::load_operation_golden("simple");
    assert_eq!(golden.turbine_powers.len(), 1);
    assert_eq!(golden.turbine_powers[0].len(), 1); // 1 turbine in this config
}

#[test]
fn test_load_timeseries_golden() {
    let golden = common::load_timeseries_golden();
    assert!(golden.timeseries.n_findex > 0);
    assert_eq!(
        golden.timeseries.wind_speeds.len(),
        golden.timeseries.n_findex
    );
}

#[test]
fn test_load_windrose_golden() {
    let golden = common::load_windrose_golden();
    assert!(golden.windrose.n_wind_directions > 0);
    assert!(golden.windrose.n_wind_speeds > 0);
}
