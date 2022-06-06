from dp_tools.core.check_model import ValidationProtocol, FlagCode


def check_green():
    code = FlagCode.GREEN
    message = "This is green"
    return {"code": code, "message": message}


def check_unhandled():
    code = FlagCode.DEV_UNHANDLED
    message = "This is developer unhandled"
    return {"code": code, "message": message}


def test_vv_components_with_explicit_run():
    """Expected to run ONLY B and B-children components"""
    vp = ValidationProtocol(run_components=["B-1"])

    with vp.component_start(
        name="A",
        description="",
    ):
        with vp.payload(payloads=[{}]):
            vp.add(check_green)

        with vp.component_start(
            name="A-1",
            description="",
        ):
            with vp.payload(payloads=[{}]):
                vp.add(check_green)

        with vp.component_start(
            name="A-2",
            description="",
        ):
            with vp.payload(payloads=[{}]):
                vp.add(check_green)

            with vp.component_start(
                name="A-2-1",
                description="",
            ):
                with vp.payload(payloads=[{}]):
                    vp.add(check_green)

    with vp.component_start(
        name="B",
        description="",
    ):
        with vp.payload(payloads=[{}]):
            vp.add(check_green)

        with vp.component_start(
            name="B-1",
            description="",
        ):
            with vp.payload(payloads=[{}]):
                vp.add(check_green)

        with vp.component_start(
            name="B-2",
            description="",
        ):
            with vp.payload(payloads=[{}]):
                vp.add(check_green)

            with vp.component_start(
                name="B-2-1",
                description="",
            ):
                with vp.payload(payloads=[{}]):
                    vp.add(check_green)

    print(vp.queued_checks(include_skipped_components=True))
    vp.run()
    assert set(vp.report(include_skipped=False)["flag_table"].index) == {
        ("ROOT", "B", "B-1")
    }


def test_vv_components_with_explicit_skip():
    """Expected to run skip B and B-children components but run A (all other components)"""

    vp = ValidationProtocol(skip_components=["B"])

    with vp.component_start(
        name="A",
        description="",
    ):
        with vp.payload(payloads=[{}]):
            vp.add(check_green)

        with vp.component_start(
            name="A-1",
            description="",
        ):
            with vp.payload(payloads=[{}]):
                vp.add(check_green)

    with vp.component_start(
        name="B",
        description="",
    ):
        with vp.payload(payloads=[{}]):
            vp.add(check_green)

        with vp.component_start(
            name="B-1",
            description="",
        ):
            with vp.payload(payloads=[{}]):
                vp.add(check_green)

        with vp.component_start(
            name="B-2",
            description="",
        ):
            with vp.payload(payloads=[{}]):
                vp.add(check_green)

            with vp.component_start(
                name="B-2-1",
                description="",
            ):
                with vp.payload(payloads=[{}]):
                    vp.add(check_green)

    print(vp.queued_checks(include_skipped_components=True))
    vp.run()
    assert set(vp.report(include_skipped=False)["flag_table"].index) == {
        ("ROOT", "A", "A-1"),
        ("ROOT", "A"),
    }


def test_vv_components_with_explicit_skip_and_run():
    """Expected to run "A" and "B-2" but skip "B-2-1" """
    vp = ValidationProtocol(run_components=["A", "B-2"], skip_components=["B-2-1"])

    with vp.component_start(
        name="A",
        description="",
    ):
        with vp.payload(payloads=[{}]):
            vp.add(check_green)

        with vp.component_start(
            name="A-1",
            description="",
        ):
            with vp.payload(payloads=[{}]):
                vp.add(check_green)

    with vp.component_start(
        name="B",
        description="",
    ):
        with vp.payload(payloads=[{}]):
            vp.add(check_green)

        with vp.component_start(
            name="B-1",
            description="",
        ):
            with vp.payload(payloads=[{}]):
                vp.add(check_green)

        with vp.component_start(
            name="B-2",
            description="",
        ):
            with vp.payload(payloads=[{}]):
                vp.add(check_green)

            with vp.component_start(
                name="B-2-1",
                description="",
            ):
                with vp.payload(payloads=[{}]):
                    vp.add(check_green)

    print(vp.queued_checks(include_skipped_components=True))
    vp.run()
    assert set(vp.report(include_skipped=False)["flag_table"].index) == {
        ("ROOT", "A", "A-1"),
        ("ROOT", "A"),
        ("ROOT", "B", "B-2"),
    }
