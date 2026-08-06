# Suppress logger:: output during the test run; individual tests that need
# to inspect log output (e.g. via a custom appender) still can, since they
# operate on the same logger and restore their own state via withr::defer.
old_log_threshold <- logger::log_threshold()
logger::log_threshold(logger::FATAL)
withr::defer(logger::log_threshold(old_log_threshold), teardown_env())
