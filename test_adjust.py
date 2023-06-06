from adjust import adjustInitialReferenceOffsets, adjustFinalReferenceOffsets

class TestInitial:
    def testOne(self):
        assert adjustInitialReferenceOffsets([None, None, 1, 2]) == [-1, 0, 1, 2]

    def testTwo(self):
        assert adjustInitialReferenceOffsets([None, None, 3, 4]) == [-1, -2, 3, 4]

    def testThree(self):
        assert adjustInitialReferenceOffsets([None, None, None, 5, 6, 7]) == [
            -2,
            -3,
            -4,
            5,
            6,
            7,
        ]


class TestFinal:
    def testOne(self):
        assert adjustFinalReferenceOffsets([3, 4, None, None]) == [3, 4, -5, -6]

    def testTwo(self):
        assert adjustFinalReferenceOffsets([5, 6, 7, None, None, None]) == [
            5,
            6,
            7,
            -8,
            -9,
            -10,
        ]
