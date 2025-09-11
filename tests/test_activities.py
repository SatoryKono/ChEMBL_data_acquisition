from library.activities import get_activities


def test_get_activities_limit() -> None:
    activities = get_activities(3)
    assert len(activities) == 3
    assert activities[0]["activity_id"] == 0
    assert activities[-1]["activity_id"] == 2
